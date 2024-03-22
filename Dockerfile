# Rshiny base image
FROM rocker/shiny
# Make a directory in the container
RUN mkdir /home/shiny-app

# Install R dependencies
RUN R -e "install.packages(c('shiny', 'shinythemes', 'ggplot2', 'shinyjs', 'png', 'dplyr', 'stringr'))"

# Copy the Shiny app code
COPY app.R /home/shiny-app/app.R

# Expose the application port
EXPOSE 3838

# Run the R Shiny app
CMD Rscript /home/shiny-app/app.R
