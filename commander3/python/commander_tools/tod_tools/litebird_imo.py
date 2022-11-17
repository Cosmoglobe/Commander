import litebird_sim as lbs


class LitebirdImo:

    def __init__(self, imo, version, instrument):
        self.imo = imo
        self.version = version
        self.url_prefix = f'/releases/{version}/satellite/{instrument}/'
        self.instrument = lbs.InstrumentInfo.from_imo(self.imo, f'{self.url_prefix}/instrument_info')

    def get_channel_dets(self, channel):
        return lbs.FreqChannelInfo.from_imo(self.imo, url=f'{self.url_prefix}{channel}/channel_info').detector_names

    def get_channel_names(self):
        return self.instrument.channel_names

    def get_detector_property(self, channel, detector, prop):
        detinfo = lbs.DetectorInfo.from_imo(self.imo, f'/releases/v1.3/satellite/LFT/{channel}/{detector}/detector_info')
        return getattr(detinfo, prop)
