"""
Weight tracking with EMA and innovation boost.

Tracks miner performance over time with exponential moving averages
and rewards breakthrough performance with temporary boosts.
"""

from typing import Dict, List, Any
import logging

logger = logging.getLogger(__name__)


class ScoreTracker:
    """Track scores with EMA and innovation boost."""

    def __init__(self, num_miners: int, alpha: float = 0.1, boost_factor: float = 0.15):
        """
        Initialize score tracker.

        Args:
            num_miners: Number of miners to track
            alpha: EMA smoothing factor
            boost_factor: Innovation boost multiplier
        """
        self.num_miners = num_miners
        self.alpha = alpha
        self.boost_factor = boost_factor

        # Initialize scores
        self.ema_scores = [0.0] * num_miners
        self.boost_states = [{'active': False, 'timestamp': 0, 'decay': 0.0} for _ in range(num_miners)]
        self.global_best = 0.0
        self.innovation_threshold = 0.05  # 5% improvement triggers boost

    def update(self, miner_id: int, new_score: float, timestamp: float = None) -> float:
        """
        Update score with EMA and check for innovation boost.

        Args:
            miner_id: Index of the miner
            new_score: New score to incorporate
            timestamp: Current timestamp for boost decay

        Returns:
            Updated effective score with boost
        """
        if miner_id >= self.num_miners:
            raise ValueError(f"Invalid miner_id {miner_id}")

        # Update EMA
        old_ema = self.ema_scores[miner_id]
        self.ema_scores[miner_id] = (1 - self.alpha) * old_ema + self.alpha * new_score

        # Check for innovation (new best score)
        if new_score > self.global_best * (1 + self.innovation_threshold):
            self.global_best = new_score
            self.boost_states[miner_id] = {
                'active': True,
                'timestamp': timestamp or 0,
                'decay': 1.0
            }
            logger.info(f"Miner {miner_id} achieved innovation boost with score {new_score:.4f}")

        # Calculate boost decay if active
        if self.boost_states[miner_id]['active'] and timestamp:
            elapsed = timestamp - self.boost_states[miner_id]['timestamp']
            half_life = 3600  # 1 hour half-life
            decay = 2 ** (-elapsed / half_life)
            self.boost_states[miner_id]['decay'] = decay

            # Deactivate boost if decayed too much
            if decay < 0.01:
                self.boost_states[miner_id]['active'] = False

        # Calculate effective score with boost
        boost = self.boost_states[miner_id]['decay'] if self.boost_states[miner_id]['active'] else 0
        effective_score = self.ema_scores[miner_id] * (1 + self.boost_factor * boost)

        return effective_score

    def get_normalized_weights(self) -> List[float]:
        """Get normalized weights for all miners."""
        import numpy as np

        # Get effective scores with boosts
        effective_scores = []
        for i in range(self.num_miners):
            boost = self.boost_states[i]['decay'] if self.boost_states[i]['active'] else 0
            score = self.ema_scores[i] * (1 + self.boost_factor * boost)
            effective_scores.append(score)

        # Normalize to sum to 1
        total = sum(effective_scores)
        if total > 0:
            weights = [s / total for s in effective_scores]
        else:
            # Equal weights if no scores yet
            weights = [1.0 / self.num_miners] * self.num_miners

        return weights

    def get_stats(self) -> Dict[str, Any]:
        """Get current statistics."""
        return {
            'ema_scores': self.ema_scores,
            'global_best': self.global_best,
            'active_boosts': sum(1 for b in self.boost_states if b['active']),
            'normalized_weights': self.get_normalized_weights()
        }
