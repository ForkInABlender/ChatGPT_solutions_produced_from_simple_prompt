# Build stage 1: Google Cloud SDK
FROM google/cloud-sdk:latest AS gcloud

# Install additional components
RUN gcloud components install kubectl

# Build stage 2: AWS CLI
FROM amazon/aws-cli:latest AS aws

# Install kubectl
RUN curl -LO "https://storage.googleapis.com/kubernetes-release/release/$(curl -s https://storage.googleapis.com/kubernetes-release/release/stable.txt)/bin/linux/amd64/kubectl"
RUN chmod +x kubectl
RUN mv kubectl /usr/local/bin

# Build stage 3: Azure CLI
FROM mcr.microsoft.com/azure-cli:latest AS az

# Install kubectl
RUN curl -LO "https://storage.googleapis.com/kubernetes-release/release/$(curl -s https://storage.googleapis.com/kubernetes-release/release/stable.txt)/bin/linux/amd64/kubectl"
RUN chmod +x kubectl
RUN mv kubectl /usr/local/bin

# Final stage: Combine all tools into one image
FROM ubuntu:latest

# Copy Google Cloud SDK, AWS CLI, and Azure CLI binaries
COPY --from=gcloud /google-cloud-sdk /google-cloud-sdk
COPY --from=aws /usr/local/bin/aws /usr/local/bin/aws
COPY --from=az /usr/local/bin/az /usr/local/bin/az

# Set PATH environment variable to include all three tools
ENV PATH="${PATH}:/google-cloud-sdk/bin:/usr/local/bin"

# Set default command
CMD ["bash"]
