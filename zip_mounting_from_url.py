# Dylan Kenneth Eliot & GPT for google

"""
This was needed for another project. This is to make it easier by uploads of zip files to only need to mount the
 zip file by URL. This makes it easier with the `get_importer/github_import.py` code to import code from anywhere and go.
  This also saves disk space when done correctly. This also means that all those packages commonly used can be
 mounted from github as compressed modules.

The way to use them together is like so:

```
from github_import import Git_import
from zipimport import zipimporter

zip_path, _ = __import__(urllib.request.urlretrieve("http://example.com/path/to/yourfile.zip"))
zipimporter = __import__("zipimport").zipimporter(zip_path)
your_module = zipimporter.load_module("your_module")

with Git_import(username="github_username_it_is_stored_by", repo="repository_your code_is_in", branch="name_of_branch", path_to_module="path/to/user_scripted_module_by_branch.py") as user_scripted_module_by_branch:
 pass
# or something you want your code to do. It can also go after the with block. However, this area would be good for making sure your modules are being imported correctly. 

from your_module import function_that_is_needed_by_user_script_to_function # just for example..


user_scripted_module_by_branch.function()

```

This not only saves disk so you can do more than store modules locally, but also allows more flexibility.

If you're having to import and check if you've installed the necessary packages everywhere you use your code,
 well, this saves space by keeping it in as few places possible and prevents continuously installing. This also
  makes it easier to track down bugs, and keep releases organized by version without over-writing someone else's
 good code.

"""
zip_path, _ = __import__(urllib.request.urlretrieve("http://example.com/path/to/yourfile.zip"))
zipimporter = __import__("zipimport").zipimporter(zip_path)
your_module = zipimporter.load_module("your_module")
