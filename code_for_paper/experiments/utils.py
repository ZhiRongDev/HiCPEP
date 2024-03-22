import numpy as np
import logging
logging.basicConfig(format='%(message)s', level=logging.INFO)

def flip_tracks(track1_np: np.ndarray, track2_np: np.ndarray):
    if len(track1_np) != len(track2_np):
        logging.info("The length of track1_np is different with track2_np")
        logging.info(f"Length of track1_np: {len(track1_np)}")
        logging.info(f"Length of track2_np: {len(track2_np)}")

    if np.corrcoef(track1_np, track2_np)[0][1] < 0:
        track2_np = -track2_np
    return track1_np, track2_np

def flip_track_gc(track_np: np.ndarray, gc_np: np.ndarray) -> np.ndarray:
    if len(track_np) != len(gc_np):
        logging.info("The length of track_np is different with gc_np")
        logging.info(f"Length of track_np: {len(track_np)}")
        logging.info(f"Length of gc_np: {len(gc_np)}")

    if np.mean(gc_np[track_np > 0]) < np.mean(gc_np[track_np < 0]):
        track_np = -track_np
    return track_np