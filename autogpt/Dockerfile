# Dylan Kenenth Eliot & gpt-4-plugins
# This docker file creates a custom image that extracts the ubuntu image with the gpt components. It is useful for if you're looking to have
#  a copilot or a system that can aide in development, even inside docker containers. de3343/kubernetes:latest can also be substituted so it
#   it can run and use kubernetes tools. 
# 
# 
# 
FROM significantgravitas/auto-gpt
COPY --from=ubuntu:latest / /
COPY --from=significantgravitas/auto-gpt / /
RUN apt-get update
RUN apt-get install apt-utils
ENTRYPOINT python -m autogpt
