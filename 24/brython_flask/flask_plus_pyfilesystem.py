from flask import Flask, render_template
import fs


app = Flask(__name__)

# Create a virtual filesystem
vfs = fs.open_fs('mem://')

# Mount the virtual filesystem using the FUSE filesystem
fuse_mount = fs.mount('fuse://', vfs)

# Expose the virtual filesystem to the Brython script
vfs_data = fuse_mount.to_json()



@app.route('/')
def index():
    return render_template('index.html')


app.run()
