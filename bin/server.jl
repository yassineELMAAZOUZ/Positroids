using Positroids

port=parse(Int,get(ENV,"PORT","10000"))
secure_cookie=lowercase(get(ENV,"POSITROIDS_SECURE_COOKIE","true")) in ("1","true","yes","on")
service=serve_positroids_web(host="0.0.0.0",port=port,secure_cookie=secure_cookie)
println("Positroids web service listening on port $port")
wait()
