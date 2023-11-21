# Dylan Kenneth Eliot & GPT-4-plugins (Alpha Edition) & Google Bard


"""

What it does is let you use kubernetes within python code.

Thus far, it will work so long as you can manage to step around string literal evaluative mangling. Beyond that, it works within the directory of kubernetes.
After that, configure as you will, as long as the go.mod and other essentials are available for kubernetes to still function.

git clone https://github.com/kubernetes/kubernetes.git
"""


import subprocess
import ctypes
import os
import tempfile

go_code = """
package main

import "C"
import (
    "context"
    "fmt"
    "os"

    metav1 "k8s.io/apimachinery/pkg/apis/meta/v1"
    "k8s.io/client-go/kubernetes"
    "k8s.io/client-go/tools/clientcmd"
)

//export ListPods
func ListPods(namespace *C.char) {
    ns := C.GoString(namespace)

    // Path to the k3s kubeconfig file
    kubeconfigPath := "/etc/rancher/k3s/k3s.yaml"

    // Check if the kubeconfig file exists
    if _, err := os.Stat(kubeconfigPath); os.IsNotExist(err) {
        fmt.Println("kubeconfig file not found at:", kubeconfigPath)
        return
    }

    // Build the config from the kubeconfig file
    config, err := clientcmd.BuildConfigFromFlags("", kubeconfigPath)
    if err != nil {
        fmt.Println("Error building kubeconfig:", err)
        return
    }

    // Create a Kubernetes client
    clientset, err := kubernetes.NewForConfig(config)
    if err != nil {
        fmt.Println("Error creating Kubernetes client:", err)
        return
    }

    // List pods in the specified namespace
    pods, err := clientset.CoreV1().Pods(ns).List(context.TODO(), metav1.ListOptions{})
    if err != nil {
        fmt.Println("Error listing pods:", err)
        return
    }

    if len(pods.Items) == 0 {
        fmt.Println("No pods found in namespace:", ns)
        return
    }

    fmt.Println("Pods in namespace", ns, ":")
    for _, pod := range pods.Items {
        fmt.Println(pod.GetName())
    }
}

func main() {}
"""

# Path to your Go workspace
go_workspace_path = "./"
go_file_path = os.path.join(go_workspace_path, "kubectl.go")
so_file_path = go_file_path.replace(".go", ".so")

# Write the Go code to the file
with open(go_file_path, "w") as go_file:
    go_file.write(go_code)

# Compile the Go code into a shared library
compilation_result = subprocess.run(["go", "build", "-o", so_file_path, "-buildmode=c-shared", go_file_path])

# Check if the compilation was successful
if compilation_result.returncode != 0:
    print("Compilation failed.")
    exit(1)

# Check if the shared library file exists
if not os.path.exists(so_file_path):
    print("Shared library file not found.")
    exit(1)

# Load the shared library in Python
try:
    kubectl = ctypes.CDLL(so_file_path)
except OSError as e:
    print(f"Error loading shared library: {e}")
    exit(1)

# Define the argument and result types for the exported Go function
kubectl.ListPods.argtypes = [ctypes.c_char_p]
kubectl.ListPods.restype = ctypes.c_void_p

# Use the function from the shared library
kubectl.ListPods(b"acorn-system")  # Replace 'default' with the desired namespace

# Clean up: remove the Go file and the compiled library
os.remove(go_file_path)
os.remove(so_file_path)

