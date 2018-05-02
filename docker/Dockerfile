FROM python:3.5
# Install missing deps
USER root
# Install unix requirements for the docker
RUN apt-get update && apt-get install -y --no-install-recommends \
        libxml2-dev     \
        libxslt1-dev    \
        libopenblas-dev \
        liblapack-dev   \
		less			\
		build-essential \
		gfortran        \
		fort77          \
		wget            \
		cmake           \
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

# Copy tnstall sundials script
COPY install_sundials.sh .
# Install sundials
RUN ./install_sundials.sh
# Export environment vairables from installation
ENV SUNDIALS_INCLUDEDIR ="${HOME}/sundials-2.7.0/include"
ENV SUNDIALS_LIBDIR ="${HOME}/sundials-2.7.0/lib"
ENV SUNDIALS_INST="${HOME}/sundials-2.7.0"
# Install python interface to sundials
RUN pip install scikits.odes

# Make the workspace folder that will link the sources
RUN mkdir /skimpy

COPY .bashrc $HOME
RUN chown "$USER" "$HOME/.bashrc"

#Finalizing installation

USER $USER
RUN mkdir ./work
WORKDIR ./work

ENV LD_LIBRARY_PATH="${HOME}/sundials-2.7.0/lib"
# Load your package in development mode on startup
ENTRYPOINT ["/bin/bash", "-c", "pip install --user -e /skimpy[ORACLE] && $0 $*"]
CMD /bin/bash