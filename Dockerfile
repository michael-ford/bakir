FROM continuumio/miniconda3

RUN apt-get update && apt-get install -y git wget && \
    echo "Cloning kir-annotator repository..." && \
    git clone https://github.com/michael-ford/kir-annotator.git /kir-annotator && \
    echo "Creating Conda environment..." && \
    conda env create -f /kir-annotator/kir-annotator-env.yml && \
    echo "Activating environment and installing kir-annotator..." && \
    . /opt/conda/etc/profile.d/conda.sh && \
    conda activate kir-annotator-env && \
    pip install /kir-annotator

ENV PATH /opt/conda/envs/kir-annotator-env/bin:$PATH

COPY entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]

