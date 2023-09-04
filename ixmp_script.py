#Dylan Kenneth Eliot & GPT-4-plugins & openai module

"""
This is a ixmp server. I don't know it's use case beyond vaguely that it could be used for testing and analyzing of data and application usage.
How is beyond me. Even the documentation was a little vague. So, probably GPT is better for asking what it can be used for.


As this is the example it wrote entirely by a fat-finger (ICMP). 
"""

import subprocess, ixmp

subprocess.run(['ixmp', 'server', 'start'])
output = subprocess.check_output(['ixmp', 'info']).decode('utf-8')
address_port = output.split(':')[1].strip()
config = {{
    'url': 'http://{address_port}',
    'backend': 'jdbc',
    'jdbc_args': {{
        'driver': 'postgresql',
        'url': 'jdbc:postgresql://localhost:5432/ixmp',
        'user': 'user',
        'password': 'password',
        'driver_args': {{
            'sslmode': 'require'
        }}
    }}
}}

mp = ixmp.Platform(**config)
