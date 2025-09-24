# Use R base image
FROM r-base

# Make a directory in the container
RUN mkdir /home/shiny-app

# Install R dependencies
RUN R -e "install.packages(c('shiny', 'shinythemes', 'ggplot2', 'shinyjs', 'png', 'dplyr', 'stringr'), Ncpus=4, repos='https://cloud.r-project.org/')"


RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


RUN R -e "install.packages(c('curl', 'openssl', 'httr'), Ncpus=4, repos='https://cloud.r-project.org/')"

# Install BiocManager
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/')"

# Install Bioconductor packages
RUN R -e "BiocManager::install('Biostrings')"

RUN R -e "install.packages(c('tidyr'))"


# Copy the Shiny app code
COPY app-prod.R /home/shiny-app


EXPOSE 3838


