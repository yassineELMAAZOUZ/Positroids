FROM julia:1.12-bookworm

WORKDIR /app
COPY Project.toml Manifest.toml ./
COPY src ./src
COPY bin ./bin
RUN julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile(); using Positroids'

ENV PORT=10000
EXPOSE 10000
CMD ["julia", "--project=/app", "--startup-file=no", "bin/server.jl"]
