#!/bin/bash

# Add Docker GPG key
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -

# Add Docker repository
sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"

# Update package list
sudo apt-get update

# Install Docker
sudo apt-get install docker-ce -y

# Start Docker service
sudo systemctl start docker

# Verify Docker installation
sudo docker run hello-world


# chatGPT wrote all that. The robot wrote it. And I tested it. It works even on a Zorin15.3 ISO booted 
# from a cd. Literally found the thing I would have searched hours to find. The bot finds it in minutes. 
