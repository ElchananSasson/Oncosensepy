from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QWidget, QListWidget, QPushButton, QHBoxLayout, QVBoxLayout, QMessageBox
from PyQt5.QtGui import QFont


class AssignValuesWindow(QWidget):
    def __init__(self, control_list, inhibitor_list, cell_name):
        """
        This function initializes an instance of the AssignValuesWindow class.

        Params:
            control_list (list): The initial values for the left list.
            inhibitor_list (list): The initial values for the right list.
            cell_name (str): The name of the cell for which the groups are being assigned.
        """
        super().__init__()

        self.result = None
        self.resize(500, 400)
        self.setWindowTitle(cell_name + " groups")

        font = QFont('sans-serif', 10)

        # Save the default values for later use in reset
        self.control_list_default = control_list
        self.inhibitor_list_default = inhibitor_list

        # Initialize the two list widgets
        self.left_list = QListWidget()
        self.right_list = QListWidget()

        # Populate the left list with all the values
        self.left_list.addItems(control_list)
        self.right_list.addItems(inhibitor_list)

        # Set the font of each list
        self.left_list.setFont(font)
        self.right_list.setFont(font)

        # Connect the itemSelectionChanged signal of each list to a method that ensures only one item is selected at
        # a time
        self.left_list.itemClicked.connect(self.on_list_item_clicked)
        self.right_list.itemClicked.connect(self.on_list_item_clicked)

        # Create buttons to move items between the lists
        self.move_right_button = QPushButton(">")
        self.move_left_button = QPushButton("<")
        self.delete_button = QPushButton("Delete")
        self.reset_button = QPushButton("Reset")
        self.apply_button = QPushButton("Apply")

        # Connect the buttons to their corresponding functions
        self.move_right_button.clicked.connect(self.move_right)
        self.move_left_button.clicked.connect(self.move_left)
        self.delete_button.clicked.connect(self.delete)
        self.reset_button.clicked.connect(self.reset)
        self.apply_button.clicked.connect(self.apply)

        # Create horizontal and vertical layouts for the buttons and list widgets
        button_layout = QVBoxLayout()
        button_layout.addWidget(self.move_right_button)
        button_layout.addWidget(self.move_left_button)
        button_layout.addWidget(self.delete_button)
        button_layout.addWidget(self.reset_button)
        button_layout.addWidget(self.apply_button)

        list_layout = QHBoxLayout()
        list_layout.addWidget(self.left_list)
        list_layout.addLayout(button_layout)
        list_layout.addWidget(self.right_list)

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

        self.move_right_button.setStyleSheet(button_style)
        self.move_left_button.setStyleSheet(button_style)
        self.delete_button.setStyleSheet(delete_style)
        self.reset_button.setStyleSheet(button_style)
        self.apply_button.setStyleSheet(button_style)

    def on_list_item_clicked(self, item):
        """
        This method is triggered when an item in either the left or right list is clicked.
        It ensures that only one item can be selected at a time in either list.

        Params:
            item (QListWidgetItem): The item that was clicked in the list.
        """

        sender_list = self.sender()  # Get the list that sent the signal
        other_list = self.left_list if sender_list is self.right_list else self.right_list

        # Clear the selection in the other list
        other_list.clearSelection()

        # Select the clicked item in the sender list
        sender_list.setCurrentItem(item)

    def move_right(self):
        """
        This method moves the selected item(s) from the left list to the right list.
        """

        selected_items = self.left_list.selectedItems()
        if not selected_items:
            return
        for item in selected_items:
            # If the item is already in the right list, don't move it
            if self.right_list.findItems(item.text(), Qt.MatchExactly):
                continue
            # Add the item to the right list and deselect any other selected items in that list
            self.right_list.addItem(item.text())
            self.right_list.clearSelection()
            # Remove the item from the left list
            self.left_list.takeItem(self.left_list.row(item))
        # Clear the selection in the left list
        self.left_list.clearSelection()

    def move_left(self):
        """
                This method moves the selected item(s) from the right list to the left list.
                """
        # Move the selected item(s) from the right list to the left list
        selected_items = self.right_list.selectedItems()
        if not selected_items:
            return
        for item in selected_items:
            # If the item is already in the left list, don't move it
            if self.left_list.findItems(item.text(), Qt.MatchExactly):
                continue
            # Add the item to the left list and deselect any other selected items in that list
            self.left_list.addItem(item.text())
            self.left_list.clearSelection()
            # Remove the item from the right list
            self.right_list.takeItem(self.right_list.row(item))
        # Clear the selection in the right list
        self.right_list.clearSelection()

    def delete(self):
        """
            This method delete the selected item from a list.
            """
        # Remove the selected item(s) from both lists
        selected_left = self.left_list.selectedItems()
        selected_right = self.right_list.selectedItems()

        if not selected_left and not selected_right:
            # If no items are selected, show an error message
            QMessageBox.warning(self, "Error", "Please select an item to delete.")
            return

        for item in selected_left:
            self.left_list.takeItem(self.left_list.row(item))
        for item in selected_right:
            self.right_list.takeItem(self.right_list.row(item))

    def reset(self):
        """
            This method reset the values to the initial values of each list.
                """
        # Reset the lists to their initial state
        self.left_list.clear()
        self.right_list.clear()
        self.left_list.addItems(self.control_list_default)
        self.right_list.addItems(self.inhibitor_list_default)

    def apply(self):
        """
          This method set the result for the selected values.
                        """
        # Get the final state of the lists and return them as two separate lists
        left_values = [self.left_list.item(i).text() for i in range(self.left_list.count())]
        right_values = [self.right_list.item(i).text() for i in range(self.right_list.count())]
        self.result = (left_values, right_values)
        self.apply_clicked = True
        self.close()

    def closeEvent(self, event):
        """
    This method handles the close event when the user clicks on the exit button.

    Params:
        event (QCloseEvent): The close event object.

    Returns:
        None
    """
        if not self.apply_clicked:
            reply = QMessageBox.question(self, 'Confirm Exit', 'Are you sure you want to exit?',
                                         QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply == QMessageBox.Yes:
                # User confirmed exit, close the window
                self.result = (self.control_list_default, self.inhibitor_list_default)
                event.accept()
            else:
                # User cancelled exit, ignore the close event
                event.ignore()
