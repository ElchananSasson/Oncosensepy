import sys

from PyQt5.QtWidgets import QWidget, QListWidget, QPushButton, QHBoxLayout, QVBoxLayout, QMessageBox
from PyQt5.QtGui import QFont


class AssignNamesValuesWindow(QWidget):
    def __init__(self, cell_line_list):
        super().__init__()

        self.result = None
        self.resize(350, 400)
        self.setWindowTitle("Cell Names:")

        font = QFont('sans-serif', 10)

        # Save the default values for later use in reset
        self.cell_names_list_default = cell_line_list

        # Initialize the two list widgets
        self.one_list = QListWidget()

        # Populate the left list with all the values
        self.one_list.addItems(cell_line_list)

        # Set the font of each list
        self.one_list.setFont(font)

        # Create buttons to move items between the lists
        self.delete_button = QPushButton("Delete")
        self.reset_button = QPushButton("Reset")
        self.apply_button = QPushButton("Apply")

        # Connect the buttons to their corresponding functions
        self.delete_button.clicked.connect(self.delete)
        self.reset_button.clicked.connect(self.reset)
        self.apply_button.clicked.connect(self.apply)

        # Create horizontal and vertical layouts for the buttons and list widgets
        button_layout = QVBoxLayout()
        button_layout.addWidget(self.delete_button)
        button_layout.addWidget(self.reset_button)
        button_layout.addWidget(self.apply_button)

        list_layout = QHBoxLayout()
        list_layout.addWidget(self.one_list)
        list_layout.addLayout(button_layout)

        # Set the main layout for the window
        self.setLayout(list_layout)

        # Flag to indicate if Apply button was clicked
        self.apply_clicked = False

        # Set styles for the buttons
        button_style = """
                    QPushButton {
                        border: 1px solid black;
                        color: black;
                        border-radius: 5px;
                        padding: 5px;
                        font-family: sans-serif

                    }
                    QPushButton:hover {
                        background-color: #2980b9;
                        cursor:pointer;
                    }
                """
        delete_style = """
                    QPushButton {
                        border: 1px solid black;
                        color: black;
                        border-radius: 5px;
                        padding: 5px;
                        font-family: sans-serif

                    }
                    QPushButton:hover {
                        background-color: #e41b1b;
                        cursor:pointer;
                    }
                """

        self.delete_button.setStyleSheet(delete_style)
        self.reset_button.setStyleSheet(button_style)
        self.apply_button.setStyleSheet(button_style)

    def delete(self):
        # Remove the selected item(s) from list
        selected_items = self.one_list.selectedItems()

        if not selected_items:
            # If no items are selected, show an error message
            QMessageBox.warning(self, "Error", "Please select an item to delete.")
            return

        for item in selected_items:
            self.one_list.takeItem(self.one_list.row(item))

    def reset(self):
        # Reset the lists to their initial state
        self.one_list.clear()
        self.one_list.addItems(self.cell_names_list_default)

    def apply(self):
        # Get the final state of the lists and return them as two separate lists
        values = [self.one_list.item(i).text() for i in range(self.one_list.count())]
        self.result = values
        self.apply_clicked = True
        self.close()

    def closeEvent(self, event):
        if not self.apply_clicked:
            reply = QMessageBox.question(self, 'Confirm Exit', 'Are you sure you want to exit?',
                                         QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply == QMessageBox.Yes:
                # User confirmed exit, close the window
                self.result = (self.cell_names_list_default)
                event.accept()
            else:
                # User cancelled exit, ignore the close event
                event.ignore()
