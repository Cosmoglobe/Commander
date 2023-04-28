import litebird_sim as lbs


class LitebirdImo:

    def __init__(self, imo, version, instrument):
        self.imo = imo
        self.version = version
        self.url_prefix = f'/releases/{version}/satellite/{instrument}/'
        self.instrument = lbs.InstrumentInfo.from_imo(self.imo, f'{self.url_prefix}/instrument_info')
        self.instrument_name = instrument

    def get_channel_dets(self, channel):
        return lbs.FreqChannelInfo.from_imo(self.imo, url=f'{self.url_prefix}{channel}/channel_info').detector_names

    def get_channel_names(self):
        return self.instrument.channel_names

    def get_detector_property(self, channel, detector, prop):
        detinfo = lbs.DetectorInfo.from_imo(self.imo, f'/releases/{self.version}/satellite/{self.instrument_name}/{channel}/{detector}/detector_info')
        return getattr(detinfo, prop)

    def get_detector_frequency(self, channel):
        return lbs.FreqChannelInfo.from_imo(self.imo, url=f'{self.url_prefix}{channel}/channel_info').bandcenter_ghz

    def get_detector_bandwidth(self, channel):
        return lbs.FreqChannelInfo.from_imo(self.imo, url=f'{self.url_prefix}{channel}/channel_info').bandwidth_ghz

    def get_detector_fwhm(self, channel):
        return lbs.FreqChannelInfo.from_imo(self.imo, url=f'{self.url_prefix}{channel}/channel_info').fwhm_arcmin
