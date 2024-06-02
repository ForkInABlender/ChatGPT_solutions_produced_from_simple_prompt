# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
The exists as a more formal bridge between javascript land and python land of execution and embedding the web in applications. 
 Or in some cases, vise versa. As seen in this example. 


"""


from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtWebChannel import QWebChannel
from PyQt5.QtCore import pyqtSlot, QObject
import sys

class Callbacks(QObject):
    @pyqtSlot()
    def buttonClicked(self):
        print("Button was clicked!")

    @pyqtSlot(str)
    def linkClicked(self, url):
        print(f"Link clicked: {url}")

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Interactive HTML Renderer with Callbacks")
        self.setGeometry(0, 0, 1024, 768)

        # Set up the central widget and layout
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        self.layout = QVBoxLayout(self.central_widget)

        # Set up the QWebEngineView
        self.browser = QWebEngineView()
        self.layout.addWidget(self.browser)

        # Define the HTML content with JavaScript
        html_content = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Sample HTML Page</title>
    <script type="text/javascript" src="https://cdn.jsdelivr.net/npm/brython@3.9.5/brython.min.js"></script>
    <script src="qrc:///qtwebchannel/qwebchannel.js"></script>
    <script>
        function setupChannel() {
            new QWebChannel(qt.webChannelTransport, function(channel) {
                window.pyCallback = channel.objects.pyCallback;
            });
        }
    </script>
</head>
<body onload="setupChannel(); brython();">
    <h1>Hello, World!</h1>
    <p>This is a sample HTML page generated using Python.</p>
    <a href="https://www.example.com" id="my_link">Click this link</a>
    <button id="my_button">Click Me</button>

    <script type="text/python">
        from browser import document, window, alert

        def show_alert(event):
            alert('Button clicked!')
            window.pyCallback.buttonClicked()

        def link_clicked(event):
            event.preventDefault()
            window.pyCallback.linkClicked(event.target.href)

        document["my_button"].bind("click", show_alert)
        document["my_link"].bind("click", link_clicked)
    </script>
</body>
</html>
        """

        # Load the HTML content
        self.browser.setHtml(html_content)

        # Set up the QWebChannel
        self.channel = QWebChannel()
        self.callbacks = Callbacks()
        self.channel.registerObject("pyCallback", self.callbacks)
        self.browser.page().setWebChannel(self.channel)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
