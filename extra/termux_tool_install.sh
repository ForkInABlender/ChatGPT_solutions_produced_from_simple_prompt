#Dylan Kenneth Eliot & GPT-4-plugins ( Beta Edition )

"""

This is to install acorn.io, kubectl, k3s, k3d, ubuntu in a proot container, and the minimum to use kubernetes on android. 

"""

#!/bin/bash







#!/bin/bash

# Define installation paths and versions
K3S_VERSION="v1.29.1+k3s2"
ACORN_VERSION="v0.10.1"
UDOCKER_URL="https://github.com/indigo-dc/udocker/releases/download/devel/udocker-1.3.0.tar.gz"
K3D_VERSION="v4.4.8"
K3D_URL="https://github.com/rancher/k3d/releases/download/${K3D_VERSION}/k3d-linux-arm64"

# Update and install necessary packages
echo "Updating and installing necessary packages..."
apt-get update && apt-get upgrade -y
apt-get install -y wget curl tar

# Install udocker
echo "Installing udocker..."
wget -O udocker.tar.gz $UDOCKER_URL
tar -xzf udocker.tar.gz
./udocker install
rm udocker.tar.gz

# Setup alias for udocker (optional)
alias udocker="./udocker"
alias docker="udocker"

# Download k3s binary
echo "Downloading k3s..."
wget -O k3s "https://github.com/k3s-io/k3s/releases/download/${K3S_VERSION}/k3s"
chmod +x k3s
mv k3s /usr/local/bin/

# Download and setup acorn
echo "Downloading acorn..."
wget "https://github.com/acorn-io/acorn/releases/download/${ACORN_VERSION}/acorn-${ACORN_VERSION}-linux-arm64.tar.gz"
tar -xzf "acorn-${ACORN_VERSION}-linux-arm64.tar.gz" -C /usr/local/bin/ acorn
chmod +x /usr/local/bin/acorn
rm "acorn-${ACORN_VERSION}-linux-arm64.tar.gz"

# Attempt to download k3d
echo "Attempting to download k3d..."
wget -O k3d $K3D_URL
chmod +x k3d
mv k3d /usr/local/bin/

echo "Installation of udocker, k3s, acorn, and an attempt for k3d within proot-distro is complete."
echo "Please note: Running k3d in this setup is experimental and may not work as expected due to Docker dependencies."
