import { useEffect, useState } from "react";

export function useWindowLocationOrigin(): string {
  const [uri, setUri] = useState("");
  useEffect(() => {
    // Bypass SSR and set genuine window.location.origin in browser
    setUri(window.location.origin);
  }, [])
  return uri;
}
