### Unused functions (yet)

function wavelenght_nm_to_omega_au(wavel_nm)
    result=45.5633525285/wavel_nm
    return result
end

function omega_au_wavelenght_nm(omega_au)
    result=45.5633525285/omega_au
    return result
end

function intensity_au_to_Wcm2(int_au)
    result=int_au*3.509444836958223e16
    return result
end

function intensity_Wcm2_to_au(int_Wcm2)
    result=int_Wcm2/3.509444836958223e16
    return result
end

function angle_rad_to_deg(angle_rad)
    result = (angle_rad *180.0)/(np.pi)
    return result
end

function angle_deg_to_rad(angle_deg)
    result = (angle_deg /180.0)*(np.pi)
    return result
end

function distance_atomic_units_to_Angstorm(au)
    result = au / 1.889725989
    return result
end

function distance_Angstorm_to_atomic_units(au)
    result = au * 1.889725989
    return result
end




### General function

function au_to_au(number)
    return number
end



### Energy functions
  
function au_to_eV(energy_au)
    result=energy_au*27.2113961
    return result
end

function eV_to_au(energy_eV)
    result=energy_eV/27.2113961
    return result
end

function meV_to_au(energy_meV)
    energy_eV = energy_meV * 1e-3
    result = energy_eV / 27.2113961 
    return result
end

function K_to_au(Kelvin)
    result_to_ev = Kelvin/11604.4251662
    result = eV_to_au(result_to_ev)
    return result
end 

### Time functions

function au_to_fs(time_au)
    result=time_au*0.024188843267176925
    return result
end

function fs_to_au(time_fs)
    result=time_fs/0.024188843267176925
    return result
end

### Laser functions

function au_to_Vm(electric_field_au)
    result = electric_field_au * 5.14220652e11
    return result
end

function Vm_to_au(electric_field_Vm)
    result = electric_field_Vm / 5.14220652e11
    return result
end

function au_to_Wcm2(int_au)
    result=int_au*3.509444836958223e16
    return result
end

function Wcm2_to_au(int_Wcm2)
    result=int_Wcm2/3.509444836958223e16
    return result
end

function nm_to_au(wavel_nm)
    result=45.5633525285/wavel_nm
    return result
end

function au_to_nm(omega_au)
    result=45.5633525285/omega_au
    return result
end