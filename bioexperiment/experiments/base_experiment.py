import datetime as dt
import pathlib
from abc import ABC, abstractmethod

root = pathlib.Path(__file__).parent.parent.parent


class Experiment(ABC):
    """
    Experiment class defines interface and deals with
    experiment folder setup.
    """

    def __init__(self, name: str):
        self.name = name
        self.experiment_date = dt.datetime.today()
        self.folder_name = f"{self.name}-{self.experiment_date.strftime('%Y-%m-%d')}"
        self.results_full_path = root / "results" / self.folder_name
        self._make_results_folder()

    def _make_results_folder(self):
        self.results_full_path.mkdir(exist_ok=True, parents=True)

    @abstractmethod
    def run(self):
        """
        Defines the steps necessary to run an experiment.
        """
