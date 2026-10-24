#
# Runtime environment ONLY.  The repository is intentionally NOT copied into
# the image.  Mount it at /workspace when running:
#
#   docker run -it --rm -v "$PWD":/workspace -w /workspace mutatedsequence:2.0.0
#
# The image contains Python 3.12 plus the optional command-line tools used for
# reference preparation and cross-checks (samtools, bcftools, tabix/bgzip).
# The Python pipeline itself has no third-party runtime dependencies.
FROM python:3.12-slim-bookworm

LABEL org.opencontainers.image.title="mutatedSequence runtime" \
      org.opencontainers.image.description="Python 3.12 + samtools/bcftools/tabix environment; source is bind-mounted at /workspace" \
      org.opencontainers.image.version="2.0.0"

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONPATH=/workspace \
    LANG=C.UTF-8

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      ca-certificates \
      curl \
      wget \
      less \
      procps \
      passwd \
      samtools \
      bcftools \
      tabix \
 && rm -rf /var/lib/apt/lists/*

# Run as a normal user by default.  On Linux, add
#   --user "$(id -u):$(id -g)"
# if your host UID is not 1000, so output files are owned by you.
RUN useradd --create-home --uid 1000 --shell /bin/bash mutant \
 && mkdir -p /workspace \
 && chown -R mutant:mutant /workspace /home/mutant

WORKDIR /workspace
USER mutant

# Interactive shell by default.  For a one-shot command:
#   docker run --rm ... mutatedsequence:2.0.0 \
#     python -m mutatedseq run ...
CMD ["bash"]
