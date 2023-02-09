FROM rocker/rstudio:4.2

RUN apt-get update \
  && apt-get -y install apt-utils \
  && apt-get -y install libxml2-dev libudunits2-dev libgdal-dev libgeos-dev libproj-dev \
  && apt-get clean

#ENV RENV_VERSION 0.15.5
#RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
#RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

#WORKDIR /opt/meta
#COPY renv.lock renv.lock

# approach one
#ENV RENV_PATHS_LIBRARY renv/library

# approach two
#RUN mkdir -p renv
#COPY .Rprofile .Rprofile
#COPY renv/activate.R renv/activate.R
#COPY renv/settings.dcf renv/settings.dcf

# RUN R -e "renv::restore()"
