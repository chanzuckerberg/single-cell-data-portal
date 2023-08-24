import React, {
  createContext,
  useCallback,
  useContext,
  useEffect,
  useState,
} from "react";

interface FullScreenContextProps {
  isFullScreen: boolean;
  enableFullScreen: () => void;
  disableFullScreen: () => void;
  screenDimensions: { width: number; height: number };
}

const FullScreenContext = createContext<FullScreenContextProps | undefined>(
  undefined
);

const FullScreenProvider: React.FC<React.PropsWithChildren> = ({
  children,
}) => {
  const [isFullScreen, setIsFullScreen] = useState(false);

  const [screenDimensions, setScreenDimensions] = useState<{
    width: number;
    height: number;
  }>({
    width: 0,
    height: 0,
  });

  const enableFullScreen = () => setIsFullScreen(true);
  const disableFullScreen = useCallback(() => setIsFullScreen(false), []);

  useEffect(() => {
    const handleResize = () => {
      if (typeof window !== "undefined") {
        setScreenDimensions({
          width: window.innerWidth,
          height: window.innerHeight,
        });
      }
    };
    window.addEventListener("resize", handleResize);
    if (isFullScreen) {
      handleResize();
    }
    return () => {
      window.removeEventListener("resize", handleResize);
    };
  }, [isFullScreen]);

  useEffect(() => {
    const handleKeyDown = (event: KeyboardEvent) => {
      if (event.key === "Escape") {
        disableFullScreen();
      }
    };

    window.addEventListener("keydown", handleKeyDown);
    return () => {
      window.removeEventListener("keydown", handleKeyDown);
    };
  }, [disableFullScreen]);

  return (
    <FullScreenContext.Provider
      value={{
        isFullScreen,
        enableFullScreen,
        disableFullScreen,
        screenDimensions,
      }}
    >
      {children}
    </FullScreenContext.Provider>
  );
};

export const useFullScreen = (): FullScreenContextProps => {
  const context = useContext(FullScreenContext);
  if (!context) {
    throw new Error("useFullScreen must be used within a FullScreenProvider");
  }
  return context;
};

export default FullScreenProvider;
