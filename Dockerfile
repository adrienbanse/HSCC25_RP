FROM julia:1.11.3

RUN mkdir -p /app/output
WORKDIR /app

COPY Project.toml ./
COPY RP.jl ./

CMD ["julia", "-e", "using Pkg; Pkg.activate(\".\"); Pkg.instantiate(); include(\"RP.jl\")"]