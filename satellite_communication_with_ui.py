import math
import tkinter as tk
from tkinter import scrolledtext, messagebox

class SatelliteCommunication:
    def __init__(self, frequency, distance, tx_power, tx_gain, rx_gain, noise_power, elevation_angle=30, rain_rate=0):
        """
        Initialize satellite communication parameters.
        :param frequency: Signal frequency in Hz
        :param distance: Distance between satellite and receiver in meters (slant range or direct distance)
        :param tx_power: Transmitter power in dBm
        :param tx_gain: Transmitter antenna gain in dBi
        :param rx_gain: Receiver antenna gain in dBi
        :param noise_power: Noise power in dBm
        :param elevation_angle: Elevation angle in degrees (default 30)
        :param rain_rate: Rain rate in mm/hr for attenuation calculation (default 0)
        """
        self.frequency = frequency
        self.distance = distance
        self.tx_power = tx_power
        self.tx_gain = tx_gain
        self.rx_gain = rx_gain
        self.noise_power = noise_power
        self.elevation_angle = elevation_angle
        self.rain_rate = rain_rate
        self.c = 3e8  # Speed of light in m/s
        self.G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
        self.M_earth = 5.972e24  # Earth's mass in kg
        self.R_earth = 6.371e6  # Earth's radius in meters

    def free_space_path_loss(self):
        """
        Calculate Free Space Path Loss (FSPL) in dB.
        FSPL(dB) = 20 * log10(d) + 20 * log10(f) + 20 * log10(4π/c)
        """
        fspl = 20 * math.log10(self.distance) + 20 * math.log10(self.frequency) + 20 * math.log10(4 * math.pi / self.c)
        return fspl

    def eirp(self):
        """
        Calculate Effective Isotropic Radiated Power (EIRP) in dBm.
        EIRP(dBm) = Pt(dBm) + Gt(dBi)
        """
        return self.tx_power + self.tx_gain

    def link_budget(self):
        """
        Calculate received power using the link budget equation.
        Pr(dBm) = EIRP(dBm) + Gr(dBi) - FSPL(dB) - Rain_Attenuation(dB)
        """
        fspl = self.free_space_path_loss()
        rain_attenuation = self.rain_attenuation()
        received_power = self.eirp() + self.rx_gain - fspl - rain_attenuation
        return received_power

    def signal_to_noise_ratio(self):
        """
        Calculate Signal-to-Noise Ratio (SNR) in dB.
        SNR(dB) = Pr(dBm) - Pn(dBm)
        """
        received_power = self.link_budget()
        snr = received_power - self.noise_power
        return snr

    def carrier_to_noise_ratio(self):
        """
        Calculate Carrier-to-Noise Ratio (C/N) in dB.
        C/N(dB) = Pr(dBm) - Pn(dBm) (same as SNR for simplicity)
        """
        return self.signal_to_noise_ratio()

    def antenna_gain(self, antenna_diameter, efficiency=0.6):
        """
        Calculate antenna gain in dBi.
        G(dBi) = 10 * log10(η * (π * D * f / c)^2)
        :param antenna_diameter: Antenna diameter in meters
        :param efficiency: Antenna efficiency (default 0.6)
        """
        gain = 10 * math.log10(efficiency * (math.pi * antenna_diameter * self.frequency / self.c) ** 2)
        return gain

    def orbital_period(self, altitude):
        """
        Calculate orbital period in seconds.
        T = 2π * sqrt((R_earth + h)^3 / (G * M_earth))
        :param altitude: Satellite altitude above Earth's surface in meters
        """
        semi_major_axis = self.R_earth + altitude
        period = 2 * math.pi * math.sqrt(semi_major_axis ** 3 / (self.G * self.M_earth))
        return period

    def doppler_shift(self, relative_velocity):
        """
        Calculate Doppler shift in Hz.
        Δf = (v / c) * f
        :param relative_velocity: Relative velocity between satellite and receiver in m/s
        """
        doppler_shift = (relative_velocity / self.c) * self.frequency
        return doppler_shift

    def slant_range(self, altitude):
        """
        Calculate slant range distance in meters.
        d = sqrt(R_earth^2 + (R_earth + h)^2 - 2 * R_earth * (R_earth + h) * cos(θ))
        where θ is the angle subtended at Earth's center.
        :param altitude: Satellite altitude above Earth's surface in meters
        """
        theta = math.radians(90 + self.elevation_angle)
        slant_range = math.sqrt(self.R_earth**2 + (self.R_earth + altitude)**2 - 
                               2 * self.R_earth * (self.R_earth + altitude) * math.cos(theta))
        return slant_range

    def rain_attenuation(self):
        """
        Calculate rain attenuation in dB (simplified ITU-R model).
        A = k * R^a * L, where L is the effective path length.
        """
        if self.rain_rate == 0:
            return 0.0
        k = 0.0188 if self.frequency < 12e9 else 0.0691
        a = 1.276 if self.frequency < 12e9 else 1.099
        effective_length = self.distance / math.sin(math.radians(self.elevation_angle))
        attenuation = k * (self.rain_rate ** a) * (effective_length / 1000)
        return attenuation

    def bit_error_rate(self, modulation="QPSK"):
        """
        Estimate Bit Error Rate (BER) based on SNR and modulation scheme.
        Simplified model for QPSK: BER ≈ 0.5 * erfc(sqrt(SNR_linear/2))
        """
        snr_db = self.signal_to_noise_ratio()
        snr_linear = 10 ** (snr_db / 10)
        if modulation == "QPSK":
            ber = 0.5 * math.erfc(math.sqrt(snr_linear / 2))
        else:
            raise ValueError("Only QPSK modulation is supported in this simplified model.")
        return ber

class SatelliteCommUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Satellite Communication Calculator")
        self.root.geometry("600x700")

        # Input fields
        self.labels = [
            "Frequency (GHz):", "Distance (km):", "Tx Power (dBm):", 
            "Tx Gain (dBi):", "Rx Gain (dBi):", "Noise Power (dBm):",
            "Elevation Angle (deg):", "Rain Rate (mm/hr):", 
            "Altitude (km):", "Relative Velocity (m/s):", "Antenna Diameter (m):"
        ]
        self.entries = {}
        self.defaults = [2.4, 1000, 30, 20, 18, -90, 30, 10, 500, 7500, 1.0]

        for i, (label, default) in enumerate(zip(self.labels, self.defaults)):
            tk.Label(root, text=label).grid(row=i, column=0, padx=5, pady=5, sticky="e")
            entry = tk.Entry(root)
            entry.insert(0, str(default))
            entry.grid(row=i, column=1, padx=5, pady=5)
            self.entries[label] = entry

        # Output text area
        self.output_text = scrolledtext.ScrolledText(root, height=15, width=60)
        self.output_text.grid(row=len(self.labels), column=0, columnspan=2, padx=5, pady=10)

        # Calculate button
        tk.Button(root, text="Calculate", command=self.calculate).grid(row=len(self.labels)+1, column=0, columnspan=2, pady=10)

    def calculate(self):
        try:
            # Get input values
            freq = float(self.entries["Frequency (GHz):"].get()) * 1e9  # Convert to Hz
            dist = float(self.entries["Distance (km):"].get()) * 1e3  # Convert to meters
            tx_power = float(self.entries["Tx Power (dBm):"].get())
            tx_gain = float(self.entries["Tx Gain (dBi):"].get())
            rx_gain = float(self.entries["Rx Gain (dBi):"].get())
            noise_power = float(self.entries["Noise Power (dBm):"].get())
            elevation_angle = float(self.entries["Elevation Angle (deg):"].get())
            rain_rate = float(self.entries["Rain Rate (mm/hr):"].get())
            altitude = float(self.entries["Altitude (km):"].get()) * 1e3  # Convert to meters
            relative_velocity = float(self.entries["Relative Velocity (m/s):"].get())
            antenna_diameter = float(self.entries["Antenna Diameter (m):"].get())

            # Initialize SatelliteCommunication
            sat_comm = SatelliteCommunication(freq, dist, tx_power, tx_gain, rx_gain, noise_power, elevation_angle, rain_rate)

            # Calculate results
            results = [
                f"EIRP: {sat_comm.eirp():.2f} dBm",
                f"Free Space Path Loss: {sat_comm.free_space_path_loss():.2f} dB",
                f"Rain Attenuation: {sat_comm.rain_attenuation():.2f} dB",
                f"Received Power: {sat_comm.link_budget():.2f} dBm",
                f"SNR: {sat_comm.signal_to_noise_ratio():.2f} dB",
                f"C/N: {sat_comm.carrier_to_noise_ratio():.2f} dB",
                f"Antenna Gain (1m dish): {sat_comm.antenna_gain(antenna_diameter):.2f} dBi",
                f"Orbital Period: {sat_comm.orbital_period(altitude) / 3600:.2f} hours",
                f"Doppler Shift: {sat_comm.doppler_shift(relative_velocity):.2f} Hz",
                f"Slant Range: {sat_comm.slant_range(altitude) / 1000:.2f} km",
                f"Bit Error Rate (QPSK): {sat_comm.bit_error_rate():.2e}"
            ]

            # Display results
            self.output_text.delete(1.0, tk.END)
            self.output_text.insert(tk.END, "\n".join(results))

        except ValueError as e:
            messagebox.showerror("Error", "Please enter valid numerical values.")
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred: {str(e)}")

def main():
    root = tk.Tk()
    app = SatelliteCommUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()