const TOAST_STORAGE_PREFIX = "byod-toast-dismissed-";

// Forced re-trigger dates
const FORCED_RETRIGGER_DATES = [
  new Date("2025-11-19"), // SCverse Workshop
  new Date("2025-12-01"), // NeurIPS (week of Dec 1, 2025)
];

// 2 months in milliseconds
const TWO_MONTHS_MS = 2 * 30 * 24 * 60 * 60 * 1000;

export interface ToastDismissalState {
  toastId: string;
  dismissedAt: number;
  dismissed: boolean;
  lastForcedDatePassed?: string; // Latest forced date that has passed (YYYY-MM-DD)
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
    if (!state.dismissed) return false; // Not dismissed, so show toast

    const now = Date.now();
    const timeSinceDismissal = now - state.dismissedAt;
    const today = new Date();

    // Check if it's been more than 2 months since dismissal
    if (timeSinceDismissal > TWO_MONTHS_MS) {
      return false; // Show toast again after 2 months
    }

    // Find the latest forced date that has passed
    const passedForcedDates = FORCED_RETRIGGER_DATES.filter(
      (date) => date <= today
    ).sort((a, b) => b.getTime() - a.getTime()); // Sort descending

    if (passedForcedDates.length > 0) {
      const latestPassedDate = passedForcedDates[0];
      const latestPassedDateString = latestPassedDate
        .toISOString()
        .split("T")[0];
      const lastForcedDatePassed = state.lastForcedDatePassed || "";

      // If a new forced date has passed since last dismissal, show toast
      if (latestPassedDateString > lastForcedDatePassed) {
        return false; // Show toast due to new forced date
      }
    }

    return true; // Keep dismissed
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
    const today = new Date();

    // Find the latest forced date that has passed
    const passedForcedDates = FORCED_RETRIGGER_DATES.filter(
      (date) => date <= today
    ).sort((a, b) => b.getTime() - a.getTime()); // Sort descending

    const latestPassedDate =
      passedForcedDates.length > 0 ? passedForcedDates[0] : null;
    const lastForcedDatePassed = latestPassedDate
      ? latestPassedDate.toISOString().split("T")[0]
      : undefined;

    const state: ToastDismissalState = {
      toastId,
      dismissedAt: Date.now(),
      dismissed: true,
      lastForcedDatePassed,
    };

    localStorage.setItem(
      `${TOAST_STORAGE_PREFIX}${toastId}`,
      JSON.stringify(state)
    );
  } catch (error) {
    console.warn("Error saving toast dismissal state:", error);
  }
}
