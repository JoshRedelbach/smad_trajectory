import math

# Constants
g0 = 9.81  # Earth's gravitational acceleration (m/s^2)

# Rocket equation to compute propellant mass for a given delta-v, I_sp, and final mass
def propellant_mass(delta_v, I_sp, mf):
    m0 = math.exp(delta_v / (I_sp * g0)) * mf
    return m0 - mf

# Inputs (example values, modify these)
delta_v1 = 7000  # Delta-v for Earth escape (m/s)
delta_v2 = 100   # Delta-v for Jupiter flyby (m/s)
delta_v3 = 4800  # Delta-v for Saturn orbit insertion (m/s)

I_sp_stage1 = 450  # Specific impulse of the first engine (s)
I_sp_stage2 = 312  # Specific impulse of the second engine (s)

m_tank_stage_1 = 1050
m_dry_stage_1 = 350 + m_tank_stage_1;        # Mass of first engine/stage (kg)

m_tank_stage_2 = 100
m_sc = 500 + m_tank_stage_2;                # mass of spacecraft kg assumption

# Third stage: Saturn orbit insertion
m_stage3_final = m_sc                 # Final mass at this stage
propellant3 = propellant_mass(delta_v3, I_sp_stage2, m_stage3_final)

# Second stage: Jupiter flyby (initial mass includes propellant for third stage)
m_stage2_final = m_sc + propellant3
propellant2 = propellant_mass(delta_v2, I_sp_stage2, m_stage2_final)

# First stage: Earth escape (initial mass includes propellant for second and third stages)
m_stage1_final = m_sc + propellant2 + propellant3 + m_dry_stage_1
propellant1 = propellant_mass(delta_v1, I_sp_stage1, m_stage1_final)

# Total mass at liftoff
total_mass_at_liftoff = m_sc + propellant1 + propellant2 + propellant3 + m_dry_stage_1
total_mass_of_propellant = propellant1 + propellant2 + propellant3

# Print results
print()
print(f"Total mass at liftoff: {total_mass_at_liftoff:.2f} kg")
print(f"Propellant mass for Earth escape: {propellant1:.2f} kg")
print(f"Propellant mass for Jupiter flyby: {propellant2:.2f} kg")
print(f"Propellant mass for Saturn orbit insertion: {propellant3:.2f} kg")
print(f"Propellant mass in total: {total_mass_of_propellant:.2f} kg")
print()


# Check if propellant would fit into tank capacity
cap_required_1 = 0.06 * propellant1
cap_required_2 = 0.04 * (propellant3 + propellant2)

if cap_required_1 <= m_tank_stage_1:
    print("Tank 1 is sufficient.")
else: 
    print("Tank 1 is not sufficient.")

if cap_required_2 <= m_tank_stage_2:
    print("Tank 2 is sufficient.")
else: 
    print("Tank 2 is not sufficient.")

print()