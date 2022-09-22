const EVENT_LIST = [
    # Event 3 in Bjorge-Engeland et al.
    (name = "event3",
     L = 462 * co.kilo, Ipeak = 466.0 * co.kilo, rise = 12.4e-6, decay = 8.2e-6),

    # Event 4 in Bjorge-Engeland et al
    (name = "event4",
     L = 84 * co.kilo, Ipeak = 257.0 * co.kilo, rise = 11.6e-6, decay = 7e-6),

    # For these events I don't have rise and decay times, so just setting
    # 12 and 8 us for all of them.
    # Ex. 1 in the notes
    (name = "ex1",
     L = 75 * co.kilo, Ipeak = 79.0 * co.kilo, rise = 12e-6, decay = 8e-6),
    
    (name = "ex2",
     L = 190 * co.kilo, Ipeak = 79.0 * co.kilo, rise = 12e-6, decay = 8e-6),

    (name = "ex3",
     L = 110 * co.kilo, Ipeak = 82.0 * co.kilo, rise = 12e-6, decay = 8e-6),

    (name = "ex4",
     L = 120 * co.kilo, Ipeak = 72.0 * co.kilo, rise = 12e-6, decay = 8e-6),    
]
