// NOTE: The same logic is duplicated in Explorer. Please keep them in sync.
// File in Explorer: client/src/components/BYODConnectors/BYODToast/utils.ts
const TOAST_STORAGE_PREFIX = "byod-discover-toast-dismissed-";

// Forced re-trigger dates
const FORCED_RETRIGGER_DATES = [
  new Date("2025-11-19").getTime(), // SCverse Workshop
  new Date("2025-12-01").getTime(), // NeurIPS (week of Dec 1, 2025),
];

// 2 months in milliseconds
const TWO_MONTHS_MS = 2 * 30 * 24 * 60 * 60 * 1000;

export interface ToastDismissalState {
  toastId: string;
  dismissedAt: number;
}

/**
 * Check if a toast has been dismissed by the user
 * Returns false (show toast) if:
 * - Never dismissed
 * - Dismissed more than 2 months ago
 * - A forced re-trigger date has passed since last dismissal
 */
export function isToastDismissed(toastId: string): boolean {
  if (typeof window === "undefined") return false;

  try {
    const stored = localStorage.getItem(`${TOAST_STORAGE_PREFIX}${toastId}`);
    if (!stored) return false; // Never dismissed, so show toast

    const state: ToastDismissalState = JSON.parse(stored);

    const now = Date.now();
    const timeSinceDismissal = now - state.dismissedAt;

    // Check if it's been more than 2 months since dismissal
    if (timeSinceDismissal > TWO_MONTHS_MS) {
      return false; // Show toast again after 2 months
    }

    // Find any forced date between dismissal and now
    const passedForcedDates = FORCED_RETRIGGER_DATES.filter(
      (date) => date <= now
    ).filter((date) => date > state.dismissedAt);
    // If there are not any, don't show the toast again. Otherwise, show it again.
    return passedForcedDates.length === 0; // Keep dismissed if no forced dates passed
  } catch (error) {
    console.warn("Error reading toast dismissal state:", error);
    return false;
  }
}

/**
 * Mark a toast as dismissed by the user
 * Records the latest forced date that has passed to prevent re-showing until next forced date
 */
export function dismissToast(toastId: string): void {
  if (typeof window === "undefined") return;

  try {
    const state: ToastDismissalState = {
      toastId,
      dismissedAt: Date.now(),
    };
    localStorage.setItem(
      `${TOAST_STORAGE_PREFIX}${toastId}`,
      JSON.stringify(state)
    );
  } catch (error) {
    console.warn("Error saving toast dismissal state:", error);
  }
}
