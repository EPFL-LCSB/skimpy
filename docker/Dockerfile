FROM python:3.5
# Install missing deps
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
        libxml2-dev     \
        libxslt1-dev    \
		less			\
    && rm -rf /var/lib/apt/lists/*

ENV USER skimpy
ENV HOME /home/$USER

RUN useradd -ms "/bin/bash" "$USER"
USER $USER
WORKDIR $HOME

USER root

# Copy python package requirements
COPY requirements.txt .

# Install python packages
RUN pip install -r requirements.txt

# Add extra src if needed
RUN mkdir /src
COPY src/ /src/

# Make the workspace folder that will link the sources
RUN mkdir /skimpy

COPY .bashrc $HOME
RUN chown "$USER" "$HOME/.bashrc"

#Finalizing installation

USER $USER
RUN mkdir ./work
WORKDIR ./work

# Load your package in development mode on startup
ENTRYPOINT ["/bin/bash", "-c", "pip install --user -e /skimpy && $0 $*"]
CMD /bin/bash