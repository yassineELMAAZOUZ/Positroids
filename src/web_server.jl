mutable struct _WebGraphSession
    current::Base.RefValue{PlabicGraph}
    has_graph::Base.RefValue{Bool}
    history::Vector{PlabicGraph}
    future::Vector{PlabicGraph}
    guard::ReentrantLock
    last_seen::Float64
end

_WebGraphSession() = _WebGraphSession(Ref(plabic_graph([1])),Ref(false),
                                      PlabicGraph[],PlabicGraph[],
                                      ReentrantLock(),time())

"""A running multi-user Positroids web service returned by `serve_positroids_web`."""
mutable struct PositroidsWebService
    server::Any
    url::String
    sessions::Dict{String,_WebGraphSession}
    closed::Base.RefValue{Bool}
end

Base.show(io::IO,service::PositroidsWebService) =
    print(io,"PositroidsWebService(\"",service.url,"\")")

function Base.close(service::PositroidsWebService)
    service.closed[]=true
    close(service.server)
    return nothing
end

function _web_push_history!(session::_WebGraphSession,G::PlabicGraph,max_history::Int)
    push!(session.history,G)
    length(session.history)>max_history && popfirst!(session.history)
end

function _web_json_response(body;status=200)
    HTTP.Response(status,["Content-Type"=>"application/json; charset=utf-8"],body)
end

function _web_error(err;status=400)
    HTTP.Response(status,["Content-Type"=>"text/plain; charset=utf-8"],sprint(showerror,err))
end

function _web_append_header!(response,key,value)
    push!(response.headers,SubString(String(key))=>SubString(String(value)))
    return response
end

function _web_parse_permutation(raw,max_n)
    permutation=_parse_interactive_permutation(raw)
    length(permutation)<=max_n ||
        throw(ArgumentError("at most $max_n boundary vertices are allowed on the public server"))
    return permutation
end

function _web_parse_drawing(raw,max_n,max_internal,max_edges)
    input=JSON.parse(String(raw))
    n=Int(get(input,"n",0))
    vertices=get(input,"vertices",Any[])
    edges=get(input,"edges",Any[])
    n<=max_n || throw(ArgumentError("at most $max_n boundary vertices are allowed"))
    length(vertices)<=max_internal ||
        throw(ArgumentError("at most $max_internal internal vertices are allowed"))
    length(edges)<=max_edges || throw(ArgumentError("at most $max_edges edges are allowed"))
    colors=Symbol[Symbol(lowercase(string(get(v,"color","")))) for v in vertices]
    positions=[(Float64(v["x"]),Float64(v["y"])) for v in vertices]
    parsed_edges=[(Int(e[1]),Int(e[2])) for e in edges]
    return plabic_graph_from_drawing(n,colors,parsed_edges,positions;reduce=true)
end

function _web_parse_drawing_animation(raw,max_n,max_internal,max_edges)
    input=JSON.parse(String(raw))
    n=Int(get(input,"n",0))
    vertices=get(input,"vertices",Any[])
    edges=get(input,"edges",Any[])
    n<=max_n || throw(ArgumentError("at most $max_n boundary vertices are allowed"))
    length(vertices)<=max_internal ||
        throw(ArgumentError("at most $max_internal internal vertices are allowed"))
    length(edges)<=max_edges || throw(ArgumentError("at most $max_edges edges are allowed"))
    colors=Symbol[Symbol(lowercase(string(get(v,"color","")))) for v in vertices]
    positions=[(Float64(v["x"]),Float64(v["y"])) for v in vertices]
    parsed_edges=[(Int(e[1]),Int(e[2])) for e in edges]
    return _drawing_reduction_animation(n,colors,parsed_edges,positions)
end

function _web_session_response(session::_WebGraphSession,request;
                               iterations,restarts,max_n,max_internal,max_edges,
                               max_history,max_boundary_cells)
    target=String(request.target)
    if request.method=="GET" && target=="/"
        return HTTP.Response(200,["Content-Type"=>"text/html; charset=utf-8"],
                             _INTERACTIVE_PLABIC_HTML)
    elseif request.method=="GET" && target=="/state"
        body=lock(session.guard) do
            state=session.has_graph[] ?
                _interactive_state_json(session.current[];iterations=iterations,restarts=restarts) :
                _interactive_blank_state_json()
            _interactive_history_json(state,!isempty(session.history),!isempty(session.future))
        end
        return _web_json_response(body)
    elseif request.method=="POST" && target=="/permutation"
        try
            body=lock(session.guard) do
                permutation=_web_parse_permutation(String(request.body),max_n)
                replacement=plabic_graph(permutation)
                session.has_graph[] && _web_push_history!(session,session.current[],max_history)
                empty!(session.future)
                session.current[]=replacement
                session.has_graph[]=true
                _interactive_history_json(
                    _interactive_state_json(replacement;iterations=iterations,restarts=restarts),
                    !isempty(session.history),false)
            end
            return _web_json_response(body)
        catch err
            return _web_error(err)
        end
    elseif request.method=="POST" && target=="/custom-graph"
        try
            body=lock(session.guard) do
                replacement,stages=_web_parse_drawing_animation(
                    request.body,max_n,max_internal,max_edges)
                session.has_graph[] && _web_push_history!(session,session.current[],max_history)
                empty!(session.future)
                session.current[]=replacement
                session.has_graph[]=true
                state=_interactive_history_json(
                    _interactive_state_json(replacement;iterations=iterations,restarts=restarts),
                    !isempty(session.history),false)
                _interactive_animation_json(stages,state)
            end
            return _web_json_response(body)
        catch err
            return _web_error(err)
        end
    elseif request.method=="POST" && startswith(target,"/move/")
        session.has_graph[] || return _web_error(ArgumentError("enter and draw a permutation first"))
        raw=target[length("/move/")+1:end]
        label=try
            isempty(raw) ? Int[] : parse.(Int,split(raw,','))
        catch
            return _web_error(ArgumentError("invalid face label"))
        end
        try
            body=lock(session.guard) do
                embedded,cycle=_square_cycle_by_label(session.current[],label;
                                                       iterations=iterations,restarts=restarts)
                moved,stages=_square_move_animation_stages(embedded,cycle)
                _web_push_history!(session,session.current[],max_history)
                empty!(session.future)
                session.current[]=moved
                state=_interactive_history_json(
                    _interactive_state_json(moved;iterations=iterations,restarts=restarts),true,false)
                _interactive_animation_json(stages,state)
            end
            return _web_json_response(body)
        catch err
            return _web_error(err)
        end
    elseif request.method=="POST" && target=="/measurement"
        session.has_graph[] || return _web_error(ArgumentError("enter and draw a permutation first"))
        try
            body=lock(session.guard) do
                _measurement_payload(session.current[],request.body;
                                     iterations=iterations,restarts=restarts)
            end
            return _web_json_response(body)
        catch err
            return _web_error(err)
        end
    elseif request.method=="POST" && target=="/facets"
        session.has_graph[] || return _web_error(ArgumentError("enter and draw a permutation first"))
        try
            body=lock(session.guard) do
                _facets_payload(session.current[];
                                iterations=iterations,restarts=restarts)
            end
            return _web_json_response(body)
        catch err
            return _web_error(err)
        end
    elseif request.method=="POST" && target=="/f-vector"
        session.has_graph[] || return _web_error(ArgumentError("enter and draw a permutation first"))
        try
            body=lock(session.guard) do
                _f_vector_payload(session.current[];max_cells=max_boundary_cells)
            end
            return _web_json_response(body)
        catch err
            return _web_error(err)
        end
    elseif request.method=="POST" && target=="/undo"
        try
            body=lock(session.guard) do
                isempty(session.history) && throw(ArgumentError("there is no previous graph"))
                push!(session.future,session.current[])
                session.current[]=pop!(session.history)
                session.has_graph[]=true
                _interactive_history_json(
                    _interactive_state_json(session.current[];iterations=iterations,restarts=restarts),
                    !isempty(session.history),true)
            end
            return _web_json_response(body)
        catch err
            return _web_error(err)
        end
    elseif request.method=="POST" && target=="/redo"
        try
            body=lock(session.guard) do
                isempty(session.future) && throw(ArgumentError("there is no next graph"))
                _web_push_history!(session,session.current[],max_history)
                session.current[]=pop!(session.future)
                session.has_graph[]=true
                _interactive_history_json(
                    _interactive_state_json(session.current[];iterations=iterations,restarts=restarts),
                    true,!isempty(session.future))
            end
            return _web_json_response(body)
        catch err
            return _web_error(err)
        end
    end
    return HTTP.Response(404,"not found")
end

function _web_cookie_id(request)
    raw=HTTP.header(request,"Cookie","")
    for part in split(raw,';')
        pair=split(strip(part),'=';limit=2)
        length(pair)==2 && pair[1]=="positroids_session" && return pair[2]
    end
    return nothing
end

function _web_security_headers!(response)
    _web_append_header!(response,"X-Content-Type-Options","nosniff")
    _web_append_header!(response,"Referrer-Policy","no-referrer")
    _web_append_header!(response,"Permissions-Policy","camera=(), microphone=(), geolocation=()")
    _web_append_header!(response,"Content-Security-Policy",
        "default-src 'self'; script-src 'self' 'unsafe-inline' https://cdn.jsdelivr.net; "*
        "style-src 'self' 'unsafe-inline'; img-src 'self' data:; connect-src 'self'; "*
        "font-src 'self' data:; frame-ancestors 'self' https://yelmaazouz.org https://www.yelmaazouz.org")
    return response
end

"""
    serve_positroids_web(; host="0.0.0.0", port=10000, ...)

Start a production-oriented, multi-user Positroids web service. Each browser is
assigned an isolated, HTTP-only cookie session. Inactive sessions expire from
memory, and public graph/request limits are enforced. The returned service can
be stopped with `close(service)`.
"""
function serve_positroids_web(;host="0.0.0.0",port=10000,
                              iterations=700,restarts=5,
                              session_ttl=3600.0,max_sessions=200,
                              max_n=30,max_internal=250,max_edges=600,
                              max_history=30,max_body_bytes=1_000_000,
                              max_boundary_cells=250_000,
                              secure_cookie=true,verbose=false)
    sessions=Dict{String,_WebGraphSession}()
    sessions_guard=ReentrantLock()
    closed=Ref(false)
    ttl=Float64(session_ttl)
    function cleanup!()
        cutoff=time()-ttl
        lock(sessions_guard) do
            for (id,session) in collect(sessions)
                session.last_seen<cutoff && delete!(sessions,id)
            end
        end
    end
    function handler(request)
        target=String(request.target)
        if request.method=="GET" && target=="/health"
            return HTTP.Response(200,["Content-Type"=>"application/json"],"{\"status\":\"ok\"}")
        end
        length(request.body)<=max_body_bytes ||
            return HTTP.Response(413,"request body is too large")
        supplied=_web_cookie_id(request)
        created=false
        id,session=lock(sessions_guard) do
            cleanup_cutoff=time()-ttl
            for (old_id,old_session) in collect(sessions)
                old_session.last_seen<cleanup_cutoff && delete!(sessions,old_id)
            end
            if !isnothing(supplied) && haskey(sessions,supplied)
                supplied,sessions[supplied]
            else
                length(sessions)<max_sessions || return (nothing,nothing)
                created=true
                new_id=string(UUIDs.uuid4())
                sessions[new_id]=_WebGraphSession()
                new_id,sessions[new_id]
            end
        end
        isnothing(id) && return HTTP.Response(503,"the server has reached its active-session limit")
        session.last_seen=time()
        response=_web_session_response(session,request;
            iterations=iterations,restarts=restarts,max_n=max_n,
            max_internal=max_internal,max_edges=max_edges,max_history=max_history,
            max_boundary_cells=max_boundary_cells)
        if created
            cookie="positroids_session=$id; Path=/; HttpOnly; SameSite=Lax; Max-Age=$(round(Int,ttl))"
            secure_cookie && (cookie*="; Secure")
            _web_append_header!(response,"Set-Cookie",cookie)
        end
        return _web_security_headers!(response)
    end
    server=HTTP.serve!(handler,String(host),Int(port);verbose=verbose)
    display_host=host=="0.0.0.0" ? "127.0.0.1" : String(host)
    service=PositroidsWebService(server,"http://$display_host:$(Int(port))/",sessions,closed)
    @async while !closed[]
        sleep(min(300.0,max(10.0,ttl/4)))
        closed[] || cleanup!()
    end
    return service
end
