# Why use mss & matplotlib for what open-cv can do?

Matplotlib and mss are smaller, simpler, lighter-weight, and produce little to know overhead. All open-cv removes is the complexity needed to adapt it to AI for
 training on available actions it can take. pybrain3 would not be able to make use of the data as effectively if I were to take that route.
Plus, the above works for any device that can run python. 

# Why http://file://? Why not use as is? And why only one file at a time?

It was simpler to work with in reference to CORS and flask setup of a server, including using pyngrok to provide the httpfs tunnel via REST api calls.
 file:// is useful for local files but not inter-server use. To save time, one can make use of the file as a REST call rehashed against the local file://
  file pointed to. I started with a file to prevent later complication. The reason it uses http.server is due to how GPT-4o alpha built from prompt engineer
 request. This can be rectified with flask later on.
