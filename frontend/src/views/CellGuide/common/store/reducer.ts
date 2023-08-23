export interface PayloadAction<Payload> {
  type: keyof typeof REDUCERS;
  payload: Payload;
}
export interface State {
  cellGuideTitle: string | null;
  cellGuideNav: JSX.Element | null;
  skinnyMode: boolean;
  mobileSearchIsOpen: boolean;
}

// (thuang): If you have derived states based on the state, use `useMemo`
// to cache the derived states instead of putting them in the state.
export const INITIAL_STATE: State = {
  cellGuideTitle: null,
  cellGuideNav: null,
  skinnyMode: false,
  mobileSearchIsOpen: false,
};

export const REDUCERS = {
  setCellGuideTitle,
  setCellGuideNav,
  setSkinnyMode,
  setMobileSearchIsOpen,
};

export function reducer(state: State, action: PayloadAction<unknown>): State {
  const { type } = action;

  const handler = REDUCERS[type];

  if (!handler) {
    throw new Error(`Unknown action type: ${type}`);
  }

  /**
   * (thuang): Figuring out the typing here is more work than its rewards
   * Ideally we'll use Redux Toolkit's `createSlice` here, but I think that's
   * too heavy for now
   */
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  return handler(state, action as any);
}

function setCellGuideTitle(
  state: State,
  action: PayloadAction<State["cellGuideTitle"]>
): State {
  return {
    ...state,
    cellGuideTitle: action.payload,
  };
}

function setCellGuideNav(
  state: State,
  action: PayloadAction<State["cellGuideNav"]>
): State {
  return {
    ...state,
    cellGuideNav: action.payload,
  };
}

function setSkinnyMode(
  state: State,
  action: PayloadAction<State["skinnyMode"]>
): State {
  return {
    ...state,
    skinnyMode: action.payload,
  };
}

function setMobileSearchIsOpen(
  state: State,
  action: PayloadAction<State["mobileSearchIsOpen"]>
): State {
  return {
    ...state,
    mobileSearchIsOpen: action.payload,
  };
}
