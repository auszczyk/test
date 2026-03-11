
# 1. DANE WEJŚCIOWE (Zmapowane z atrybutów HydroIQ)
row_data = {
    'ID_IIP': 'PL.RDLP.1',          # Unikalny identyfikator
    'DLUGOSC': 1200.0,              # Długość rowu [m]
    'SZEROKOSC_DNA': 0.8,           # Szerokość dna [m]
    'NACH_SKARP': 1.5,              # Stosunek 1:n (np. dla 1:1.5 -> 1.5)
    'ID_GLEBOKOSC': 1.0,            # Głębokość (np. środek przedziału 1-1.5m)
    'Z_POCZATEK': 32.59,            # Rzędna dna na początku [m n.p.m.]
    'Z_KONIEC': 31.88,              # Rzędna dna na końcu [m n.p.m.]
    'OFFSET_ZASTAWKI': 1,           # Wysokość progu spiętrzenia. Przyjęta wysokość zastawki [m]
    'MANNING_N': 0.035              # Szorstkość (np. rów porośnięty trawą)
}


def generate_swmm_inp(data, filename="model_rowu.inp"):
    depth = data['ID_GLEBOKOSC']

    if depth < 1.0:
        h_zastawki = 1.0
    elif 1.0 <= depth <= 1.5:
        h_zastawki = 1.5
    else:
        h_zastawki = 1.7

    max_depth = h_zastawki

    # Budowa struktury sekcyjnej pliku .inp
    sections = {
        "[TITLE]": f"Model rowu HydroIQ z zastawka - ID: {data['ID_IIP']}",

        "[JUNCTIONS]": (
            # Węzeł początkowy i węzeł przed zastawką [1]
            f"J_1      {data['Z_POCZATEK']}  {max_depth}  0  0  0\n"
            f"J_Gate   {data['Z_KONIEC']}    {max_depth}  0  0  0"
        ),

        "[OUTFALLS]": (
            # Ujście (Outfall) poniżej rzędnej dna, by umożliwić swobodny odpływ
            f"Out_1    {data['Z_KONIEC'] - 0.1}  FREE  NO"
        ),

        "[CONDUITS]": (
            # Rów melioracyjny (trapezowy) łączący początek z węzłem zastawki
            f"C_1  J_1  J_Gate  {data['DLUGOSC']}  {data['MANNING_N']}  0  0"
        ),

        "[ORIFICES]": (
            # Definicja zastawki (Orifice) jako linku sterowalnego [2]
            # ID_ZASTAWKI  Węzeł_In  Węzeł_Out  Typ   Offset  Wsp_Wypływu
            f"Z_1  J_Gate  Out_1  SIDE  {data['OFFSET_ZASTAWKI']}  0.8"
        ),

        "[XSECTIONS]": (
            # Przekrój rowu (Trapezowy) i otworu zastawki (Prostokątny zamknięty) [6, 7]
            f"C_1  TRAPEZOIDAL  {max_depth}  {data['SZEROKOSC_DNA']}  {data['NACH_SKARP']}  {data['NACH_SKARP']}\n"
            f"Z_1  RECT_CLOSED  0.5  {data['SZEROKOSC_DNA']}  0  0"
        ),
        "[REPORT]": "INPUT NO\nCONTROLS NO\nSUBCATCHMENTS ALL\nNODES ALL\nLINKS ALL",
        "[OPTIONS]":
                    "FLOW_UNITS           CMS\n"
                    "INFILTRATION         HORTON\n"
                    "FLOW_ROUTING         DYNWAVE\n"
                    "START_DATE           03/10/2026\n"
                    "START_TIME           00:00:00\n"
                    "END_DATE             03/12/2026\n"
                    "END_TIME             12:00:00"

    }

    # Zapis do pliku tekstowego [5]
    with open(filename, 'w', encoding='utf-8') as f:
        for section_name, content in sections.items():
            f.write(f"{section_name}\n")
            f.write(f"{content}\n\n")

    print(f"Plik {filename} został wygenerowany poprawnie z definicją zastawki Z_1.")


# Uruchomienie generatora
generate_swmm_inp(row_data)