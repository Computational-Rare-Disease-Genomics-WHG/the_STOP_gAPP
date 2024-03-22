# Use R base image
FROM r-base

# Make a directory in the container
RUN mkdir /home/shiny-app

# Install R dependencies
RUN R -e "install.packages(c('shiny', 'shinythemes', 'ggplot2', 'shinyjs', 'png', 'dplyr', 'stringr'), Ncpus=4, repos='https://cloud.r-project.org/')"

# Copy the Shiny app code
COPY app.R /home/shiny-app

# Expose the application port
EXPOSE 3838