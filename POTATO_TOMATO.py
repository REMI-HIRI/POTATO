
from matplotlib.figure import Figure


def plot_TOMATO(FD):
    figure1 = Figure(figsize=(9, 7), dpi=70)
    subplot1 = figure1.add_subplot(111)

    F = FD[:, 0]
    PD_nm = FD[:, 1]

    subplot1.set_xlabel("Distance [nm]")
    subplot1.set_ylabel("Force [pN]")
    subplot1.plot(PD_nm, F, color='gray')
    subplot1.tick_params('both', direction='in')
    subplot1.set_ylim([min(F), max(F)])
    subplot1.set_xlim([min(PD_nm) - 10, max(PD_nm) + 10])

    return figure1
