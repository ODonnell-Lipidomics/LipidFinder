from builtins import object
import time


class Timer(object):
    """Object to keep track of the time since it's instantiation or previous mark to the
    current mark. A mark is a point where the time is recorded.

    Attributes:
        current (float): The time of the current mark
        indStart (float): The time of the previous mark
        start (float): The time when the Timer object was instantiated
    """

    def __init__(self):
        self.start = time.time()
        self.indStart = self.start
        self.current = self.start

    def total(self):
        """The total time since intantiation to the current mark.

        Returns:
            float: The number of seconds from instantiation to the current mark in seconds to 3 dp
        """
        return round((self.current - self.start), 3)

    def individual(self):
        """The total time from previous mark to current mark.

        Returns:
            float: The number of seconds from the previous mark to the current mark in seconds to 3 dp
        """
        return round((self.current - self.indStart), 3)

    def mark(self):
        """Create a new time point mark
        """
        self.indStart = self.current
        self.current = time.time()
