import pickle
import os

class Adapt:
    def __init__(self):
        self.data = [1, 2, 3]  # Daten direkt in der Initialisierung setzen

    def save_data(self, filename, folder="stored_adapt-vqe_runs"):
        """
        Speichert Daten mit pickle. Falls die Datei existiert, wird sie überschrieben.
        Der Ordner muss bereits existieren.
        
        Arguments:
            filename (str): Der Name der Datei.
            folder (str): Das Verzeichnis, in dem die Datei gespeichert wird.
        """
        if not os.path.exists(folder):  # Prüfen, ob der Ordner existiert
            raise FileNotFoundError(f"Ordner '{folder}' existiert nicht. Bitte erst erstellen.")

        filepath = os.path.join(folder, filename)  # Pfad zur Datei erstellen
        with open(filepath, "wb") as f:
            pickle.dump(self.data, f)

        print(f"Daten wurden in '{filepath}' gespeichert.")

# Beispiel:
adapt_instance = Adapt()  # Instanz der Klasse erstellen
adapt_instance.save_data("data.pkl")  # Daten speichern


