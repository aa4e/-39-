from saleae.analyzers import HighLevelAnalyzer, AnalyzerFrame, ChoicesSetting

MASK37 = 0x1FFFFFFFFF  # 37 бит


class Hla(HighLevelAnalyzer):
    """
    BiSS-C для ИДМ39А-4096-MT-BISS. Вход: SPI-анализатор (MISO = линия SLO).

    Кадр 37 бит (нумерация по PDF, MSB = 35):
      [36]=1  [35]=Ack(0)  [34]=Start(1)  [33]=CDS
      [32..20]=TURNS(12)  [19..8]=ANGLE(12)  [7]=E(инв.)  [6]=W  [5..0]=CRC6
    CRC: полином X^6+X^1+X^0 (0x43), инвертированный, считается по битам 32..6.
    """

    alignment = ChoicesSetting(choices=('Auto (поиск 1-0-1 + CRC)', 'Fixed (отбросить 4 бита)'))
    crc_check = ChoicesSetting(choices=('On', 'Off'))

    result_types = {
        'biss_frame': {
            'format': '{{data.status}} Turns={{data.turns}}; Angle={{data.angle_deg}}° {{data.flags}}'
        },
        'biss_bad': {
            'format': 'BiSS: {{data.reason}}'
        },
    }

    # Пауза между посылками. По PDF таймаут t_out = 1..20 мкс,
    # значит разрыв больше 1ms — это точно новый кадр.
    GAP_S = 1e-3 #1 ms

    def __init__(self):
        self._bits = []
        self._start = None
        self._last_end = None

    def _reset(self):
        self._bits = []
        self._start = None

    def decode(self, frame: AnalyzerFrame):
        if frame.type == 'enable':
            self._reset()
            return None

        if frame.type != 'result':
            return None

        raw_bytes = frame.data.get('miso')
        if raw_bytes is None:
            return None
        if isinstance(raw_bytes, int):
            raw_bytes = bytes([raw_bytes])

        for byte in raw_bytes:
            for i in range(7, -1, -1):
                self._bits.append((byte >> i) & 1)
                if self._start is None:
                    self._start = frame.start_time

        # Когда накопили достаточно битов, пытаемся найти валидный кадр
        if len(self._bits) >= 40:
            result = self._process(frame.end_time)
            if result is not None:
                # Нашли валидный кадр — сбрасываем буфер
                self._reset()
                return result
            else:
                # Не нашли — сдвигаемся на 1 бит и пробуем снова
                self._bits.pop(0)
                if len(self._bits) == 0:
                    self._start = None

        return None

    # ---------- разбор кадра ----------

    def _process(self, end_time):
        raw40 = 0
        for b in self._bits[-40:]:
            raw40 = (raw40 << 1) | b

        candidates = []
        if self.alignment.startswith('Fixed'):
            candidates = [raw40 & MASK37]
        else:
            for shift in (4, 3, 2, 1, 0):
                candidates.append((raw40 >> shift) & MASK37)

        chosen = None
        for cand in candidates:
            prefix_ok = (((cand >> 36) & 1) == 1 and
                         ((cand >> 35) & 1) == 0 and
                         ((cand >> 34) & 1) == 1)
            crc_ok = ((cand & 0x3F) == self._crc6(cand))
            if prefix_ok and (crc_ok or self.crc_check == 'Off'):
                chosen = cand
                break
        
        if chosen is None:
            # Кадр невалиден — возвращаем None, чтобы сдвинуться на бит
            return None

        raw = chosen
        msb    = (raw >> 36) & 1
        ack    = (raw >> 35) & 1
        start  = (raw >> 34) & 1
        cds    = (raw >> 33) & 1
        turns  = (raw >> 20) & 0x1FFF
        angle  = (raw >> 8) & 0xFFF
        e_bit  = (raw >> 7) & 1
        w_bit  = (raw >> 6) & 1
        crc_rx = raw & 0x3F
        crc_calc = self._crc6(raw)

        prefix_ok = (msb == 1 and ack == 0 and start == 1)
        crc_ok = (crc_rx == crc_calc)
        has_error = (e_bit == 0)
        has_warning = (w_bit == 1)
        angle_deg = angle * 360.0 / 4096.0

        flags = ''
        if has_error:
            flags += '[E]'
        if has_warning:
            flags += '[W]'
        #if self.crc_check == 'On' and not crc_ok:
            flags += '[CRC rx={:02X} calc={:02X}]'.format(crc_rx, crc_calc)

        ok = prefix_ok and (crc_ok or self.crc_check == 'Off')

        return AnalyzerFrame('biss_frame', self._start, end_time, {
            'status': 'OK' if ok else 'ERR',
            'turns': turns,            
            'angle_deg': round(angle_deg, 2),
            'flags': flags,
            'raw_hex': '{:09X}'.format(raw),
            'raw_bin': '{:037b}'.format(raw),
            'cds': cds,
            'e_bit': e_bit,
            'w_bit': w_bit,
            'crc_rx': '{:02X}'.format(crc_rx),
            'crc_calc': '{:02X}'.format(crc_calc),
            'crc_ok': 'yes' if crc_ok else 'no',
            'prefix_ok': 'yes' if prefix_ok else 'no',
        })

    @staticmethod
    def _crc6(raw37: int) -> int:
        POLY = 0x03
        MASK = 0x3F 
        crc = 0
        for i in range(32, 5, -1):
            bit = (raw37 >> i) & 1
            feedback = ((crc >> 5) & 1) ^ bit
            crc = (crc << 1) & MASK
            if feedback:
                crc ^= POLY

        return crc ^ 0x3F
