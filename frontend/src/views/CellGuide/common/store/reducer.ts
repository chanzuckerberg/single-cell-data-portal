export interface PayloadAction<Payload> {
  type: keyof typeof REDUCERS;
  payload: Payload;
}
export interface State {
  cellGuideTitle: string | null;
  cellGuideNav: JSX.Element | null;
}

// (thuang): If you have derived states based on the state, use `useMemo`
// to cache the derived states instead of putting them in the state.
export const INITIAL_STATE: State = {
  cellGuideTitle: null,
  cellGuideNav: null,
};

export const REDUCERS = {
  setCellGuideTitle,
  setCellGuideNav,
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
  console.log(action.payload, "******************");
  const newState = {
    ...state,
    cellGuideTitle: action.payload,
  };
  console.log(newState);
  return newState;
}

function setCellGuideNav(
  state: State,
  action: PayloadAction<State["cellGuideNav"]>
): State {
  console.log(action.payload, "******************");
  const newState = {
    ...state,
    cellGuideNav: action.payload,
  };
  console.log(newState);
  return newState;
}
