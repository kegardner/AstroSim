from dataclasses import dataclass
import scopesim as sim

@dataclass
class Command:
    sed="sne/sn1a"
    instrument ="MICADO"
    amp = 18
    filterCurve = "Ks"
    coordinates = (0,0)
    dims = (12288,12288)
    fieldOfView = "IMG_4mas"
    adaptiveOptics = "SCAO"
    airmass = 1.2 #kg/m3?
    fw1 = "open" #TODO add full name
    fw2 = "Ks" #TODO likewise
    humidity = 0.1
    temperature = 7 #C?
    pressure = 0.755 #atm?
    exposureTime = 1000 #s?

    def __repr__(self):
        return f"amp:{self.amp}"
    
    def asCMD(self):
        sim.UserCommands(use_instrument=self.instrument,
                    set_modes=[self.adaptiveOptics,self.fieldOfView],  
                    properties={"!OBS.dit": self.exposureTime,
                                "!OBS.airmass": self.airmass,
                                "!OBS.filter_name_fw1": self.fw1,
                                "!OBS.filter_name_fw2": self.fw2,
                                "!ATMO.humidity" : self.humidity,
                                "!ATMO.temperature" : self.temperature,
                                "!ATMO.pressure" : self.pressure,
                                "!DET.height" : self.dims[0],
                                "!DET.width" : self.dims[1]
                                })

    def _with(self, sed = sed, instrument = instrument, amp = amp, filterCurve = filterCurve, coordinates = coordinates, dims = dims, fieldOfView = fieldOfView, adaptiveOptics = adaptiveOptics, airmass = airmass, fw1 = fw1, fw2 = fw2, humidity = humidity, temperature = temperature, pressure = pressure, exposureTime = exposureTime):
        cmd = Command()
        cmd.sed = sed
        cmd.instrument = instrument
        cmd.amp = amp
        cmd.filterCurve = filterCurve
        cmd.coordinates = coordinates
        cmd.dims = dims
        cmd.fieldOfView = fieldOfView
        cmd.adaptiveOptics = adaptiveOptics
        cmd.airmass = airmass
        cmd.fw1 = fw1
        cmd.fw2 = fw2
        cmd.humidity = humidity
        cmd.temperature = temperature
        cmd.pressure = pressure
        cmd.exposureTime = exposureTime
        return cmd

if __name__ == "__main__":
    t = Command()
    t.amp = 25
    g = Command()
    print("Basic reflection of a python dataclass (or most classes, really...):")
    for property in filter(lambda prop : not "__" in prop,t.__dir__()):
        if t.__getattribute__(property) != g.__getattribute__(property):
            print(f"{property}: something fun")
        
        #print(f"{property}: {t.__getattribute__(property)}")