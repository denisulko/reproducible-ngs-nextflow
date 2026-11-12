FROM mambaorg/micromamba:1.5.8

LABEL maintainer="Denis Ulko <https://github.com/denisulko>" \
      description="Containerized Nextflow DSL2 Germline NGS Alignment, Variant Calling & MultiQC Telemetry Environment" \
      version="1.0.0"

USER root
WORKDIR /workspace

COPY environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

ENV PATH="/opt/conda/bin:${PATH}"
COPY . /workspace

CMD ["python", "scripts/run_multiqc_aggregator.py"]
