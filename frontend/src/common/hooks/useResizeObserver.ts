// Core dependencies
import { RefObject, useCallback, useEffect, useRef, useState } from "react";

export type ElementRect = {
  bottom: number;
  height: number;
  left: number;
  right: number;
  scrollHeight: number;
  top: number;
  width: number;
  x: number;
  y: number;
};

/**
 * Element resizing and repositioning observer.
 * @param ref - element to be observed for changes to its size or position.
 * @returns Element size and position properties for the given element.
 */
export function useResizeObserver(
  ref: RefObject<HTMLElement>
): ElementRect | undefined {
  const [elementRect, setElementRect] = useState<ElementRect>();
  const observerRef = useRef<ResizeObserver>();
  const observer = observerRef && observerRef.current;
  const observedEl = ref && ref.current;

  // Observed element is resized or repositioned - set state elementRect with the element's new dimensions or position.
  const onResize = useCallback((entries: ResizeObserverEntry[]) => {
    if (entries && entries.length > 0) {
      const entry = entries[0]; // grab the first entry; observing a single element
      const scrollHeight = entry.target.scrollHeight;
      const contentRect = entry.contentRect.toJSON();
      setElementRect({ ...contentRect, scrollHeight });
    }
  }, []);

  // Creates a new ResizeObserver object which can be used to report changes to an element's dimensions or position.
  useEffect(() => {
    observerRef.current = new ResizeObserver(onResize);
  }, [onResize]);

  // Element is "observed" and reports any changes to the element's dimensions or position.
  useEffect(() => {
    if (!observer || !observedEl) return;
    observer.observe(observedEl);
    return () => {
      observer.unobserve(observedEl);
    };
  }, [observedEl, observer]);

  return elementRect;
}
