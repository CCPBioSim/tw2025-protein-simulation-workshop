# Start with BioSim base image.
ARG BASE_IMAGE=latest
FROM ghcr.io/ccpbiosim/jupyterhub-base:$BASE_IMAGE

LABEL maintainer="James Gebbie-Rayet <james.gebbie@stfc.ac.uk>"

ARG AMBERTOOLS_DL=null

USER root
RUN mkdir /opt/amber22 && \
    chown 1000:100 -R /opt/amber22 

WORKDIR /tmp

RUN wget $AMBERTOOLS_DL/Amber22.tar.bz2 && \
    wget $AMBERTOOLS_DL/AmberTools22.tar.bz2

RUN tar xvjf AmberTools22.tar.bz2 && \
    tar xvjf Amber22.tar.bz2 && \
    rm AmberTools22.tar.bz2 Amber22.tar.bz2

WORKDIR /tmp/amber22_src
RUN ./update_amber --update
RUN mkdir /tmp/build
WORKDIR /tmp/build

RUN pip install matplotlib scipy
RUN cmake /tmp/amber22_src -DCMAKE_INSTALL_PREFIX=/opt/amber22 -DBUILD_PYTHON=TRUE -DDISABLE_TOOLS="pytraj" -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE -DDOWNLOAD_MINICONDA=FALSE -DCOMPILER=MANUAL -DBUILD_GUI=FALSE -DCOMPILER=GNU -DOPENMP=TRUE -DCUDA=FALSE -DMPI=FALSE -DUSE_FFT=True -DBUILD_DEPRECATED=False -DBUILD_INDEV=False -DBUILD_PERL=True -DOPTIMIZE=True
RUN make -j8
RUN make install

# Cleanup.
RUN rm -r /tmp/amber22_src /tmp/build

# Switch to jovyan user.
USER $NB_USER
WORKDIR $HOME

# Add the exec to the path.
ENV AMBERHOME=/opt/amber22
ENV PATH="$PATH:/opt/amber${AMBER_VERSION}/bin"
ENV LD_LIBRARY_PATH="$AMBERHOME/lib"
ENV PERL5LIB="$AMBERHOME/lib/perl"

# Install workshop deps
RUN mamba install -c conda-forge openmm ambertools
RUN pip install git+https://github.com/CompBioAsia/CBAtools.git

# Get workshop files and move them to jovyan directory.
COPY --chown=1000:100 . .
RUN rm -rf INSTALL.txt LICENSE README.md docker .git .github .gitignore

# Copy updated lab workspace
#COPY --chown=1000:100 docker/default-37a8.jupyterlab-workspace /home/jovyan/.jupyter/lab/workspaces/default-37a8.jupyterlab-workspace

# UNCOMMENT THIS LINE FOR REMOTE DEPLOYMENT
COPY docker/jupyter_notebook_config.py /etc/jupyter/

# Always finish with non-root user as a precaution.
USER $NB_USER
