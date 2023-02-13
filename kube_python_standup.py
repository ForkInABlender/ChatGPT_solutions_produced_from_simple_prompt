"""
Basic script for running and hosting your website for free. Basically is only using kubernetes to do the top-level-domain name listener and should respond by returning
 your docker container by the instance that is mapped to the IP.
"""

from kubernetes import client, config

# Load the kubeconfig file
config.load_kube_config(config_file='~/.kube/config')

# Connect to the Kubernetes API
v1 = client.CoreV1Api()

# Create a namespace
namespace = client.V1Namespace(
    api_version='v1',
    kind='Namespace',
    metadata=client.V1ObjectMeta(
        name='example-namespace'
    )
)
v1.create_namespace(body=namespace)

# Create a deployment
deployment = client.V1Deployment(
    api_version='apps/v1',
    kind='Deployment',
    metadata=client.V1ObjectMeta(
        name='example-deployment',
        namespace='example-namespace'
    ),
    spec=client.V1DeploymentSpec(
        replicas=3,
        selector={
            'matchLabels': {
                'app': 'example-app'
            }
        },
        template=client.V1PodTemplateSpec(
            metadata=client.V1ObjectMeta(
                labels={
                    'app': 'example-app'
                }
            ),
            spec=client.V1PodSpec(
                containers=[
                    client.V1Container(
                        name='example-container',
                        image='nginx:latest',
                        ports=[
                            client.V1ContainerPort(container_port=80)
                        ]
                    )
                ]
            )
        )
    )
)
apps_v1 = client.AppsV1Api()
apps_v1.create_namespaced_deployment(
    namespace='example-namespace',
    body=deployment
)

# Create a service
service = client.V1Service(
    api_version='v1',
    kind='Service',
    metadata=client.V1ObjectMeta(
        name='example-service',
        namespace='example-namespace'
    ),
    spec=client.V1ServiceSpec(
        selector={
            'app': 'example-app'
        },
        ports=[
            client.V1ServicePort(
                name='http',
                port=80,
                target_port=80
            )
        ],
        type='LoadBalancer',
        external_name='example.com'
    )
)
v1.create_namespaced_service(
    namespace='example-namespace',
    body=service
)

# Print the status of the deployment
deployment_status = apps_v1.read_namespaced_deployment_status(
    name='example-deployment',
    namespace='example-namespace'
)
print(deployment_status)
