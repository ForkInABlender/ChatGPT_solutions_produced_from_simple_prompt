# Dylan Kenneth Eliot

echo """

This script is used to reload the main server from disk.

""" >> /dev/null

apt-get update && apt-get install -y python3-pip python3.10
npm install -g localtunnel
python3.10 -m pip install --force-reinstall ./install_wheels/*.whl # flask==3.1.2 
#python3.10 -m pip install --force-reinstall flask==3.1.2
#python3.10 -m pip install -U setuptools "gradient<3.0"
python3.10 -m pip install --force-reinstall ./install_wheels/*.tar.gz



cp scipy-correcteds/scipy/__init__.py /usr/local/lib/python3.10/dist-packages/scipy/
#python3.10 AI_neurons.py
#python3.10 reload.py


run_app(){
    python3.10 app.py;
}

run_tunnel(){
    lt --port 5000 --subdomain='awakened-coffers' --print-requests &
}

run_tunnel; run_app;
