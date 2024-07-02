import { Notification } from "src/views/Collection/hooks/useNotification/types";

export const EDIT_DATASET_ERROR_NOTIFICATION: Omit<Notification, "id"> = {
  intent: "negative",
  slideDirection: "left",
  title: "Error encountered when renaming dataset.",
};
