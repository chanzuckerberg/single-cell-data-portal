const TOAST_STORAGE_PREFIX = "bottom-toast-dismissed-";

export interface ToastDismissalState {
  toastId: string;
  dismissedAt: number;
  dismissed: boolean;
}

/**
 * Check if a toast has been dismissed by the user
 */
export function isToastDismissed(toastId: string): boolean {
  if (typeof window === "undefined") return false;

  try {
    const stored = localStorage.getItem(`${TOAST_STORAGE_PREFIX}${toastId}`);
    if (!stored) return false;

    const state: ToastDismissalState = JSON.parse(stored);
    return state.dismissed;
  } catch (error) {
    console.warn("Error reading toast dismissal state:", error);
    return false;
  }
}

/**
 * Mark a toast as dismissed by the user
 */
export function dismissToast(toastId: string): void {
  if (typeof window === "undefined") return;

  try {
    const state: ToastDismissalState = {
      toastId,
      dismissedAt: Date.now(),
      dismissed: true,
    };

    localStorage.setItem(
      `${TOAST_STORAGE_PREFIX}${toastId}`,
      JSON.stringify(state)
    );
  } catch (error) {
    console.warn("Error saving toast dismissal state:", error);
  }
}

/**
 * Clear dismissal state for a toast (useful for testing or re-enabling toasts)
 */
export function clearToastDismissal(toastId: string): void {
  if (typeof window === "undefined") return;

  try {
    localStorage.removeItem(`${TOAST_STORAGE_PREFIX}${toastId}`);
  } catch (error) {
    console.warn("Error clearing toast dismissal state:", error);
  }
}
