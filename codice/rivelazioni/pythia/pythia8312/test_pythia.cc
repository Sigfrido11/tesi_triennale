#include "Pythia8/Pythia.h"

int main() {
    Pythia8::Pythia pythia;  // Crea un oggetto Pythia
    pythia.readString("Init:showProcesses = on");  // Imposta la configurazione
    pythia.init();  // Inizializza il generatore
    pythia.stat();  // Stampa le statistiche
    return 0;
}

