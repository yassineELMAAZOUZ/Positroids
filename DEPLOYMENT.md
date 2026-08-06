# Deploying the Positroids web service

The public application is a multi-user Julia HTTP service. Browser sessions are
isolated with an HTTP-only cookie and expire from server memory after one hour.

## Run locally

```bash
POSITROIDS_SECURE_COOKIE=false PORT=10000 julia --project=. --startup-file=no bin/server.jl
```

Visit `http://127.0.0.1:10000`. On Render, omit
`POSITROIDS_SECURE_COOKIE`; secure cookies are enabled by default. You can also
start and stop a local service from Julia directly:

```julia
using Positroids
service = serve_positroids_web(
    host="127.0.0.1",
    port=10000,
    secure_cookie=false,
)
```

Close a locally started service with `close(service)`.

## Render

1. Push this directory to a Git repository.
2. In Render, create a **Web Service** from that repository.
3. Choose the **Docker** runtime and `./Dockerfile`.
4. Set the health-check path to `/health`.
5. Deploy. `render.yaml` contains the equivalent service configuration.
6. Test the generated `onrender.com` URL in two private browser windows; each
   must receive an independent blank workspace.
7. Add `positroids.yelmaazouz.org` under the service's custom domains.
8. In Squarespace DNS, add the CNAME target displayed by Render using the host
   `positroids`. Do not modify the root or `www` records used by Netlify.

## Default public limits

- 200 active sessions
- one-hour inactivity expiration
- 30 boundary vertices
- 250 internal vertices
- 600 edges
- 30 undo states per session
- 1 MB request bodies

These defaults can be changed through keywords to `serve_positroids_web`.
