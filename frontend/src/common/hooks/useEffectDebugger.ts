// https://stackoverflow.com/a/59843241
import { EffectCallback, useEffect, useRef } from "react";

const usePrevious = (value: unknown, initialValue: unknown): unknown[] => {
  const ref = useRef(initialValue);
  useEffect(() => {
    ref.current = value;
  });
  return ref.current as unknown[];
};
const useEffectDebugger = (
  effectHook: EffectCallback,
  dependencies: unknown[],
  dependencyNames = []
) => {
  const previousDeps = usePrevious(dependencies, []);

  const changedDeps = dependencies.reduce(
    (
      acc: Record<string, { before: unknown; after: unknown }>,
      dependency,
      index
    ) => {
      if (dependency !== previousDeps[index]) {
        const keyName = dependencyNames[index] || index;
        return {
          ...acc,
          [keyName]: {
            before: previousDeps[index],
            after: dependency,
          },
        };
      }

      return acc;
    },
    {} as Record<string, { before: unknown; after: unknown }>
  );

  if (Object.keys(changedDeps).length) {
    console.log("[use-effect-debugger] ", changedDeps);
  }

  useEffect(effectHook, [effectHook, ...dependencies]);
};
export default useEffectDebugger;
