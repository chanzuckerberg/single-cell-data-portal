import React, { createContext, useContext, useState, ReactNode } from "react";

interface BYODModalContextType {
  isOpen: boolean;
  openModal: () => void;
  closeModal: () => void;
}

const BYODModalContext = createContext<BYODModalContextType | undefined>(
  undefined
);

interface BYODModalProviderProps {
  children: ReactNode;
}

export const BYODModalProvider: React.FC<BYODModalProviderProps> = ({
  children,
}) => {
  const [isOpen, setIsOpen] = useState(false);

  const openModal = () => setIsOpen(true);
  const closeModal = () => setIsOpen(false);

  return (
    <BYODModalContext.Provider value={{ isOpen, openModal, closeModal }}>
      {children}
    </BYODModalContext.Provider>
  );
};

export const useBYODModal = (): BYODModalContextType => {
  const context = useContext(BYODModalContext);
  if (context === undefined) {
    throw new Error("useBYODModal must be used within a BYODModalProvider");
  }
  return context;
};
