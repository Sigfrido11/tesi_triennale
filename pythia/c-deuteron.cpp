#include "Pythia8/Pythia.h"
#include <iostream>
#include <fstream>

using namespace Pythia8;

int main() {
    // Configurazione di Pythia
    Pythia pythia;

    // Configura per collisioni pesanti (Pb-Pb)
    pythia.readString("HeavyIon:mode = 1"); // Modalità collisioni ioniche
    pythia.readString("Beams:idA = 208");   // Nucleo Pb (A=208)
    pythia.readString("Beams:idB = 208");   // Nucleo Pb (A=208)
    pythia.readString("Beams:eA = 2760");   // Energia del fascio A (GeV per nucleone)
    pythia.readString("Beams:eB = 2760");   // Energia del fascio B (GeV per nucleone)
    pythia.readString("Beams:frameType = 1"); // Sistema del laboratorio

    // Aggiungi la produzione di particelle specifiche
    pythia.readString("SoftQCD:all = on"); // Attiva produzione soft QCD (dominante in Pb-Pb)

    // Inizializza Pythia
    pythia.init();

    // File di output per registrare i vertici
    std::ofstream outfile("vertex_production_c_deuterons.txt");

    // Contatore per i c-deuteroni
    int cDeuteronCount = 0;

    // Loop sugli eventi
    const int nEvents = 10000; // Numero di eventi da simulare
    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
        if (!pythia.next()) continue;

        // Analisi dell'evento
        for (int i = 0; i < pythia.event.size(); ++i) {
            // Seleziona c-deuteroni (id PDG=1000010020)
            if (pythia.event[i].id() == 1000010020) {
                ++cDeuteronCount;

                // Recupera il vertice di produzione
                double x = pythia.event[i].xProd();
                double y = pythia.event[i].yProd();
                double z = pythia.event[i].zProd();
                double t = pythia.event[i].tProd();

                // Salva nel file di output
                outfile << "Event: " << iEvent
                        << " Particle: " << pythia.event[i].name()
                        << " Vertex: (" << x << ", " << y << ", " << z << ", " << t << ")"
                        << std::endl;
            }
        }
    }

    // Statistiche finali
    pythia.stat();

    // Chiudi il file di output
    outfile.close();

    std::cout << "Simulation complete. c-Deuterons found: " << cDeuteronCount << std::endl;
    return 0;
}

