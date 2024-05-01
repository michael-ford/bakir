FROM continuumio/miniconda3

# Set the working directory in the container
WORKDIR /app

# Install git
RUN apt-get update && apt-get install -y git

# Clone the repository
RUN git clone https://github.com/michael-ford/bakir.git

# Create the conda environment
RUN conda env create -f /app/bakir/bakir-env.yml

# Initialize conda in bash shell
RUN conda init bash

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "bakir-env", "/bin/bash", "-c"]

# Install bakir within the conda environment
RUN pip install /app/bakir

# Ensure commands run inside the conda environment
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "bakir-env", "bakir"]
