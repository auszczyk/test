import numpy as np
import pandas as pd
from pyswmm import Simulation, Nodes, Links
import matplotlib.pyplot as plt

class ModelScsCn:

    def __init__(self, p_zast):
        # --- 1. PARAMETRY WEJŚCIOWE (Mapowanie z dokumentu Koncepcja HydroIQ) ---
        self.CN = 50                                # Wyznaczone z KL_PRZEP_GRUNTU i GL_UZYTK_GRUNTU
        self.A_km2 = 10.0                           # Powierzchnia zlewni (Area) - Obszar oddziaływania
        self.L_m = 1200                             # Długość rowu (DLUGOSC)
        self.I_pct = 0.06                           # Spadek dna rowu (SPADEK_DNA) [%]
        self.dt_h = 0.166                           # Krok czasowy hietogramu (15 min = 0.25h)
        self.lambda_param = 0.2                     # Parametr strat początkowych [1]
        self.k_reservoir = 0.2                      # Współczynnik opróżniania dla zlewni rolniczej
        # self.P_series = np.array([50.0])          # Przykładowy hietogram opadu P jednostkowego [mm]
        self.P_series = self.random_rain()          # Przykładowy hietogram opadu P [mm]
        self.model_rowu = r'model_rowu.inp'         # plik modelu .inp
        self.result = self.base_calculation(p_zast, self.P_series)
        self.param1 = p_zast

    def random_rain(self):
        # Parametry czasowe
        dni = 2
        krok_min = 15
        dt_h = krok_min / 60.0                      # 0.25h
        liczba_krokow = int((dni * 24) / dt_h)      # 192 kroki

        # Tworzymy tablicę zer (brak deszczu)
        self.P_series = np.zeros(liczba_krokow)

        rng = np.random.default_rng(seed=42)

        # Zakres 20:32 (12 elementów) - rozkład gamma
        zakres1 = 20
        zakres2 = 32
        self.P_series[zakres1:zakres2] = rng.gamma(shape=1.0, scale=1.1, size=zakres2-zakres1)
        print(f'max Opadu I:  {round(max(self.P_series[zakres1:zakres2]),2)}')
        print(f'min Opadu I:  {round(min(self.P_series[zakres1:zakres2]),2)}')

        zakres1 = 80
        zakres2 = 88
        self.P_series[zakres1:zakres2] = rng.gamma(shape=1.0, scale=1.1, size=zakres2-zakres1)
        print(f'max Opadu II:  {round(max(self.P_series[zakres1:zakres2]),2)}')
        print(f'min Opadu II:  {round(min(self.P_series[zakres1:zakres2]),2)}')


        zakres1 = 150
        zakres2 = 158
        self.P_series[zakres1:zakres2] = rng.uniform(0.1, 3.0, size=zakres2-zakres1)
        print(f'max Opadu III:  {round(max(self.P_series[zakres1:zakres2]),2)}')
        print(f'min Opadu III:  {round(min(self.P_series[zakres1:zakres2]),2)}')


        self.P_series = np.round(self.P_series, 2)

        return self.P_series

    # def rain_visual(self, param):

    # --- 2. OBLICZENIA HYDROLOGICZNE (Logika) ---
    def calculate_effective_rainfall(self, P_series, CN, lam):
        S_mm = (25400 / CN) - 254
        Ia_max = lam * S_mm
        print('Ia_max', Ia_max)
        print('S_mm', S_mm)

        P_cum = np.cumsum(P_series)
        Ia_dynamic = np.minimum(Ia_max, P_cum)              # Poprawka dynamiczna Ia

        # Skumulowany opad netto
        P_eff_cum = np.where(P_cum > Ia_dynamic,
                             ((P_cum - Ia_dynamic) ** 2) / (P_cum - Ia_dynamic + S_mm),
                             0)
        # print('P_eff_cum', P_eff_cum)

        # Przyrostowy opad netto
        P_eff_inc = np.diff(P_eff_cum, prepend=0)
        return P_eff_inc, S_mm


    def get_unit_hydrograph(self, A_km2, L_m, I_pct, S_mm, dt_h):
        # Konwersja jednostek dla wzoru Lag Time
        L_feet = L_m * 3.28084
        S_inches = S_mm / 25.4

        # Lag time tl [h]
        tl = (L_feet ** 0.8 * (S_inches + 1) ** 0.7) / (1900 * np.sqrt(I_pct))
        Tp = (dt_h / 2) + tl
        qp = (0.208 * A_km2) / Tp       # Przepływ szczytowy [m3/s]
        Tb = 2.67 * Tp

        # Tworzenie rzędnych trójkątnego UH
        time_steps = np.arange(0, Tb + dt_h, dt_h)
        uh = np.where(time_steps <= Tp,
                      qp * (time_steps / Tp),
                      qp * (Tb - time_steps) / (Tb - Tp))
        uh = np.maximum(uh, 0)

        # Renormalizacja objętości (Walidacja)
        V_expected = A_km2 * 1000                       # Objętość dla 1mm opadu [m3]
        V_uh = np.sum(uh) * dt_h * 3600                 # Objętość pod wykresem [m3]

        return uh * (V_expected / V_uh)


    def base_calculation(self, p_zast, P_series):
        # Obliczenia bazowe
        P_eff, S_mm = self.calculate_effective_rainfall(self.P_series, self.CN, self.lambda_param)
        uh = self.get_unit_hydrograph(self.A_km2, self.L_m, self.I_pct, S_mm, self.dt_h)

        # Splot (Convolution) - Odpływ bezpośredni Qd
        Qd = np.convolve(P_eff, uh)

        # Odpływ bazowy (Baseflow - Zbiornik liniowy)
        Qb = np.zeros(len(self.P_series))

        Infiltration = self.P_series - P_eff        # I(t)
        for t in range(1, len(self.P_series)):
            Qb[t] = Qb[t - 1] * np.exp(-self.k_reservoir * self.dt_h) + Infiltration[t] * (1 - np.exp(-self.k_reservoir * self.dt_h))

        # Rozszerzenie Qb zerami do długości Qd
        Qb_extended = np.zeros(len(Qd))
        Qb_extended[:len(Qb)] = Qb
        Q_total = Qd + Qb_extended

        # print('Q_total', Q_total)

        # --- 3. INTEGRACJA Z PYSWMM (Symulacja hydrauliczna) ---
        results = []
        with Simulation(self.model_rowu) as sim:
            node_inlet = Nodes(sim)["J_1"]                  # Węzeł wlotowy rowu
            sim.step_advance(300)                           # Krok symulacji 5 min

            zastawka = Links(sim)["Z_1"]                    # Identyfikator zastawki z pliku .inp

            for step in sim:
                aktualny_poziom = node_inlet.depth
                if p_zast in (0,1):
                    zastawka.target_setting = p_zast
                else:
                    if aktualny_poziom > 1:
                        zastawka.target_setting = 1
                        # print((round(aktualny_poziom,1), 1))
                    elif aktualny_poziom > 0.75:
                        zastawka.target_setting = 0.75
                        # print((round(aktualny_poziom,1), 0.75))
                    elif aktualny_poziom > 0.5:
                        zastawka.target_setting = 0.25
                        # print((round(aktualny_poziom,1), 0.5))
                    else:
                        zastawka.target_setting = 0
                        # print((round(aktualny_poziom,1), 1))

                # Obliczenie przepływu dla aktualnego czasu symulacji

                # Oblicz całkowity czas symulacji w godzinach
                elapsed_hours = (sim.current_time - sim.start_time).total_seconds() / 3600.0
                current_idx = int(elapsed_hours / self.dt_h)


                # Wstrzykuj dopływ tylko wtedy, gdy indeks mieści się w tablicy Q_total
                if 0 <= current_idx < len(Q_total):
                    inflow_value = Q_total[current_idx]
                else:
                    # Jeśli tablica się skończyła (indeks >= 20), woda przestała dopływać
                    inflow_value = 0.0

                # Implemenctacja (dynamiczna) przepływu do modelu
                node_inlet.generated_inflow(inflow_value)

                # Odczyt wyników w czasie rzeczywistym
                results.append({
                    'time': sim.current_time,
                    'depth': node_inlet.depth,                  # Poziom wody w rowie
                    'inflow': inflow_value
                })

            return results


if __name__ == '__main__':

    results1 = ModelScsCn(p_zast=1).result
    results2 = ModelScsCn(p_zast=0).result
    results3 = ModelScsCn(p_zast=2).result

    # --- 4. WIZUALIZACJA WYNIKÓW (Matplotlib) ---
    df_res = pd.DataFrame(results1)
    df_res2 = pd.DataFrame(results2)
    df_res3 = pd.DataFrame(results3)

    plt.figure(figsize=(10, 6))
    plt.plot(df_res['time'], df_res['depth'], label='zastawka otwarta')
    plt.plot(df_res2['time'], df_res2['depth'], label='zastawka zamknięta')
    plt.plot(df_res3['time'], df_res3['depth'], label='zastawka zmienna')

    plt.title("Reakcja rowu na odpływ wyliczony metodą hybrydową SCS-CN i SWMM\n")

    plt.xlabel("Czas")
    plt.ylabel("Głębokość [m]")
    plt.legend(title=f"Poziom wody w rowie [m]")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

