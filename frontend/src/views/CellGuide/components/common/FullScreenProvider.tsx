import React, {
  createContext,
  useCallback,
  useContext,
  useEffect,
  useState,
} from "react";
import { RIGHT_SIDEBAR_WIDTH_PX } from "../CellGuideCard/constants";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";

interface FullScreenContextProps {
  isFullScreen: boolean;
  enableFullScreen: () => void;
  disableFullScreen: () => void;
  screenDimensions: { width: number; height: number };
}

const FullScreenContext = createContext<FullScreenContextProps | undefined>(
  undefined
);

interface FullScreenProviderProps extends React.PropsWithChildren {
  cellInfoSideBarDisplayed?: boolean;
}

const FullScreenProvider: React.FC<FullScreenProviderProps> = ({
  children,
  cellInfoSideBarDisplayed,
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
          width:
            window.innerWidth -
            Number(cellInfoSideBarDisplayed) * RIGHT_SIDEBAR_WIDTH_PX,
          height: window.innerHeight - HEADER_HEIGHT_PX,
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
  }, [isFullScreen, cellInfoSideBarDisplayed]);

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
